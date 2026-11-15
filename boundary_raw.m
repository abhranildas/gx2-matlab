function [raw_grad,raw_hess,at_cusp,cusp_scale,inversion_grad_hess,nx]=boundary_raw(x,mu,v,quad,wrt,want_hess,AbsTol,RelTol)

% BOUNDARY_RAW Shared machinery behind CDF_GRAD_BD and
% NORM_ERR_GRAD_BD: one class's raw (pre-cusp-resolution, pre-squeeze,
% trailing-nx) gradient/Hessian blocks for F(x0)=P(q(x)<=x0), plus
% everything needed to resolve a cusp there -- the at_cusp mask, its scale,
% and a function handle inversion_grad_hess(s_eff,want_hess) giving the
% grad/Hessian with the normal term replaced by an arbitrary s_eff.
% CDF_GRAD_BD resolves the cusp on this class alone; NORM_ERR_GRAD_BD
% instead combines two classes' raw output and inversion probes *before*
% resolving, since resolving each class's infinity independently and then
% subtracting can throw away a genuine inf-inf cancellation in the combined
% error.
%
% Internal helper, not meant to be called directly -- see
% CDF_GRAD_BD for the documented public interface.
%
% Abhranil Das
% Center for Perceptual Systems, University of Texas at Austin
% Comments, questions, bugs to abhranil.das@utexas.edu

mu=mu(:);
q1c=quad.q1(:);
x=x(:)';                 % row of thresholds
nx=numel(x);
d=numel(mu);
wanted=@(g) any(strcmpi(g,wrt));
want_q12=wanted('q1') || wanted('q2');
opts={'AbsTol',AbsTol,'RelTol',RelTol};

% convert to gx2 params, and reuse the standardized-quadratic eigen-structure
% (S=Sigma^{1/2}, V, and the full eigenvalues d of S*Q2*S) for the per-node
% M^{-1}(t)=S*V*diag(1/(1-2i*t*d_j))*V'*S -- no d-by-d inverse per node.
[w,k,l,s,m,aux]=norm_quad_to_gx2_params(mu,v,quad);
S=aux.S; V=aux.V; dvals=aux.d(:).';     % dvals is 1-by-d
SigInv_mu=v\mu;                          % Sigma^{-1}*mu (real d-vector)

% At s=0 with mixed-sign weights, the density has a cusp at x0=m (e.g. a
% logarithmic divergence for two k=1 terms of opposite sign), and some
% derivative blocks are genuinely infinite there. A boundary built from
% exact symmetry (e.g. an antipodal-mean, swapped-covariance classification
% boundary) lands exactly on it, not just in principle. Detected once here;
% resolved by the caller via resolve_cusp_terms, which compares the
% inversion at a shrinking sequence of tiny artificial normal terms rather
% than trusting the s=0 integral right at that point.
CUSP_REL_TOL=1e-8;                    % deliberately tiny -- see gx2_derivatives.md
mixed_sign=~(all(w>0)||all(w<0));
cusp_scale=max(1,sum(abs(w).*k));
at_cusp=(s==0) & mixed_sign & (abs(x-m)<=CUSP_REL_TOL*cusp_scale);

inversion_grad_hess=@inversion_grad_hess_fn;
[g,raw_hess]=inversion_grad_hess_fn(s,want_hess);
raw_grad=struct();
for key={'q0','q1','q2'}
    if wanted(key{1}), raw_grad.(key{1})=g.(key{1}); end
end

    function [gr,hr]=inversion_grad_hess_fn(s_eff,want_hess_probe)
        % Every block is one Fourier-inversion integral of the form
        % Re(amplitude*exp(-i t (x0-m))), the amplitude being built from the
        % weights M^{-1}(t) and mu_tilde(t) times the characteristic
        % function (with its offset m moved into the kernel), integrated by
        % gx2_osc_integral. The gradient blocks are density-type integrals
        % (the (it) of the master formula cancels the 1/t of the inversion):
        % dF/dq1 = -(1/pi)\int Re[mu_tilde phi e^{-i t x0}] dt and
        % dF/dQ2 = -(1/pi)\int Re[(M^{-1}+mu_tilde mu_tilde') phi e^{-i t x0}] dt.
        % The Hessian blocks are (1/pi)\int t*Im[W(t) phi e^{-i t x0}] dt
        % for block weights W built from M^{-1} and mu_tilde, and
        % t*Im(z)=Re(-i*t*z). Parametrized by s_eff so the same code serves
        % the actual s and the cusp probe (s_eff -> tiny artificial values,
        % s itself being 0 there).
        [~,var]=gx2stat(w,k,l,s_eff,m);
        sd=sqrt(var); if sd==0, sd=1; end
        c=(x-m)/sd;                     % frequencies in the variable v=sd*t
        % dF/dq0 = -f(x0), since q0 shifts q rigidly
        gr=struct('q0',reshape(-gx2pdf(x,w,k,l,s_eff,m,opts{:}),1,nx));
        if want_q12
            A=-gx2_osc_integral(@(vv) amp(@grad_amp,vv,sd,s_eff),c,AbsTol,RelTol)/pi;
            A=reshape(A,d,d+1,nx);
            gr.q1=reshape(A(:,1,:),d,nx);
            gr.q2=symm_pages(reshape(A(:,2:d+1,:),d,d,nx));
        end
        if ~want_hess_probe, hr=[]; return; end
        Hraw=gx2_osc_integral(@(vv) amp(@hess_amp,vv,sd,s_eff),c,AbsTol,RelTol)/pi;
        off=0;
        q0q0=Hraw(1,:);                              off=1;
        q0q1=Hraw(off+(1:d),:);                      off=off+d;
        q0q2=reshape(Hraw(off+(1:d^2),:),d,d,nx);    off=off+d^2;
        q1q1=reshape(Hraw(off+(1:d^2),:),d,d,nx);    off=off+d^2;
        q1q2=reshape(Hraw(off+(1:d^3),:),d,d,d,nx);  off=off+d^3;
        q2q2=reshape(Hraw(off+(1:d^4),:),d,d,d,d,nx);
        hr=struct('q0q0',q0q0,'q0q1',q0q1,'q0q2',symm_pages(q0q2),...
            'q1q1',symm_pages(q1q1),'q1q2',q1q2,'q2q2',q2q2);
    end

    function out=amp(f,vv,sd,s_eff)
        % an amplitude f(t,s_eff) (a column) at each entry of the row vv of
        % the dimensionless variable v=sd*t, as gx2_osc_integral expects
        out=f(vv(1)/sd,s_eff)/sd;
        out=[out, zeros(numel(out),numel(vv)-1)];
        for j=2:numel(vv), out(:,j)=f(vv(j)/sd,s_eff)/sd; end
    end

    function out=grad_amp(t,s_eff)
        [Minv,mut,phi]=weights(t,s_eff);
        out=reshape([mut, Minv+mut*mut.']*phi,[],1);  % [q1-weight | Q2-weight]
    end

    function out=hess_amp(t,s_eff)
        [Minv,mut,phi]=weights(t,s_eff);
        P=Minv+mut*mut.';                 % M^{-1}+mu_tilde*mu_tilde'
        Pmm=mut*mut.';
        Wq1q2=zeros(d,d,d);
        for a=1:d, for b=1:d, for cc=1:d
            Wq1q2(a,b,cc)=2*Minv(a,b)*mut(cc)+mut(a)*P(b,cc);
        end, end, end
        WQ2=zeros(d,d,d,d);
        for a1=1:d, for a2=1:d, for b1=1:d, for b2=1:d
            WQ2(a1,a2,b1,b2)=2*(Minv(a2,b1)*Minv(b2,a1)+Minv(a2,b1)*Pmm(b2,a1) ...
                +Pmm(a2,b1)*Minv(b2,a1))+P(a2,a1)*P(b2,b1);
        end, end, end, end
        out=-1i*t*phi*[1; mut(:); P(:); P(:); Wq1q2(:); WQ2(:)];   % Nt-by-1
    end

    function [Minv,mut,phi]=weights(t,s_eff)
        % per-node characteristic function (offset m moved into the kernel)
        % and the tilted covariance/mean, from the shared eigen-structure:
        % M^{-1}=S*V*diag(1/(1-2i t d_j))*V'*S.
        phi=gx2char(t,w,k,l,s_eff,0);
        gj=1./(1-2i*t*dvals);             % 1-by-d
        Minv=S*(V.*gj)*V.'*S;
        mut=Minv*(SigInv_mu+1i*t*q1c);    % mu_tilde
    end

    function A3=symm_pages(A3)
        % symmetrize each d-by-d page (kills round-off in the plain matrix blocks)
        for pp=1:size(A3,3), A3(:,:,pp)=0.5*(A3(:,:,pp)+A3(:,:,pp).'); end
    end
end
