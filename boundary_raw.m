function [raw_grad,raw_hess,at_cusp,cusp_scale,inversion_grad_hess,nx]=boundary_raw(x,mu,v,quad,wrt,want_hess,AbsTol,RelTol,precision,n_ruben)

% BOUNDARY_RAW Shared machinery behind CDF_GRAD_NORM_QUAD and
% NORM_ERR_GRAD_BD: one class's raw (pre-cusp-resolution, pre-squeeze,
% trailing-nx) gradient/Hessian blocks for F(x0)=P(q(x)<=x0), plus
% everything needed to resolve a cusp there -- the at_cusp mask, its scale,
% and a function handle inversion_grad_hess(s_eff,want_hess) giving the
% s~=0 inversion route's grad/Hessian at an arbitrary artificial normal
% term. CDF_GRAD_NORM_QUAD resolves the cusp on this class alone;
% NORM_ERR_GRAD_BD instead combines two classes' raw output and inversion
% probes *before* resolving, since resolving each class's infinity
% independently and then subtracting can throw away a genuine inf-inf
% cancellation in the combined error.
%
% Internal helper, not meant to be called directly -- see
% CDF_GRAD_NORM_QUAD for the documented public interface.
%
% Abhranil Das
% Center for Perceptual Systems, University of Texas at Austin
% Comments, questions, bugs to abhranil.das@utexas.edu

mu=mu(:);
q1c=quad.q1(:);
x=x(:)';                 % row of thresholds
nx=numel(x);
d=numel(mu);
Nt=1+d+d^2+d^2+d^3+d^4;  % q0q0, q0q1, q0q2, q1q1, q1q2, q2q2
wanted=@(g) any(strcmpi(g,wrt));
opts={'AbsTol',AbsTol,'RelTol',RelTol,'precision',precision};

% convert to gx2 params, and reuse the standardized-quadratic eigen-structure
% (S=Sigma^{1/2}, V, and the full eigenvalues d of S*Q2*S) for the per-node
% M^{-1}(t)=S*V*diag(1/(1-2i*t*d_j))*V'*S -- no d-by-d inverse per node.
[w,k,l,s,m,aux]=norm_quad_to_gx2_params(mu,v,quad);
S=aux.S; V=aux.V; dvals=aux.d(:).';     % dvals is 1-by-d
SigInv_mu=v\mu;                          % Sigma^{-1}*mu (real d-vector)

% s=0 is the pure-quadratic-boundary (classification) regime, where the
% inversion integrals lose convergence for total dof D<=4. There we take the
% robust shifted-dof route: expand the weights in the
% eigenbasis, M^{-1}=sum_j g_j u_j u_j' and mu_tilde=sum_j g_j c_j u_j (with
% u_j=Sigma^{1/2} v_j, g_j=(1-2i t w_j)^{-1}, c_j=alpha_j+i t beta_j), so that
% every block collapses to a finite sum of shifted-dof density derivatives
% f^(n)_[..], evaluated robustly by gx2_dens_deriv. A single engine (evalblock)
% serves the gradient (one factor it, p0=1) and the Hessian (two, p0=2). When
% s~=0 (or s=0 with large D) the Gaussian damping makes the direct inversion
% converge, and we keep it.
s0=(s==0);

% The mixed-sign, s=0 route above carries x-derivatives on whichever split
% part is sampled in its smooth interior, away from its own edge (see
% gx2_dens_deriv's comments). That shielding breaks down when the threshold
% x0 coincides exactly with m, the alignment point of the two floors -- a
% boundary built from exact symmetry (e.g. an antipodal-mean, swapped-
% covariance classification boundary) lands exactly here, not just in
% principle (gx2_derivatives.md open item 3.7(a)). Detected once here;
% resolved by the caller via resolve_cusp_terms, which compares the
% (always-valid, cusp-independent) s~=0 inversion at a shrinking sequence of
% tiny artificial normal terms, rather than trusting whatever the s=0
% route's quad() returns right at that point.
CUSP_REL_TOL=1e-8;                    % deliberately tiny -- see gx2_derivatives.md open item 3.3
mixed_sign=~(all(w>0)||all(w<0));
cusp_scale=max(1,sum(abs(w).*k));
at_cusp=s0 & mixed_sign & (abs(x-m)<=CUSP_REL_TOL*cusp_scale);

if s0
    U=S*V;                    % columns u_j = Sigma^{1/2} v_j
    alph=U.'*SigInv_mu;       % d-by-1: alpha_j = u_j'*Sigma^{-1}*mu
    bet =U.'*q1c;             % d-by-1: beta_j  = u_j'*q1
    tol0=1e-9*max(1,max(abs(dvals)));
    compj=zeros(1,d);         % mode -> merged (w,k) component; 0 for a zero mode
    for j=1:d
        if abs(dvals(j))>tol0, [~,compj(j)]=min(abs(w-dvals(j))); end
    end
    memo=containers.Map('KeyType','char','ValueType','any');
    densopts={'AbsTol',AbsTol,'RelTol',RelTol,'n_ruben',n_ruben};
end

% q0 block: dF/dq0 = -f(x0). Since q0 shifts q rigidly, this is just -pdf.
raw_grad=struct();
if wanted('q0')
    if s0, raw_grad.q0=Dterm([],1); else, raw_grad.q0=reshape(-gx2pdf(x,w,k,l,s,m,opts{:}),1,nx); end
end

% q1 and Q2 blocks. Inversion route (s~=0): one integration over t returns
% both, since they share the weights M^{-1}(t) and mu_tilde(t); the (it) of the
% master formula cancels the 1/t of the inversion, leaving density-type
% integrals dF/dq1 = -(1/pi)\int Re[mu_tilde phi e^{-i t x0}] dt and
% dF/dQ2 = -(1/pi)\int Re[(M^{-1}+mu_tilde mu_tilde') phi e^{-i t x0}] dt.
% Robust route (s=0): the same two blocks from the eigenbasis engine at p0=1.
if wanted('q1') || wanted('q2')
    if s0
        if wanted('q1')
            raw_grad.q1=reshape(evalblock(mk(1,1,true),1,1),d,nx);
        end
        if wanted('q2')
            raw_grad.q2=symm_pages(reshape(evalblock([mk(1,[1 1],false),mk(1,[1 2],[true true])],2,1),d,d,nx));
        end
    else
        if strcmpi(precision,'basic')
            A=integral(@(t) integrand(t,s),0,inf,'ArrayValued',true,'AbsTol',AbsTol,'RelTol',RelTol);
        else
            A=vpa_integrate();
        end
        A=-A/pi;                              % [d, d+1, nx]
        if wanted('q1')
            raw_grad.q1=reshape(A(:,1,:),d,nx);
        end
        if wanted('q2')
            raw_grad.q2=symm_pages(reshape(A(:,2:d+1,:),d,d,nx));
        end
    end
end

% ---- Hessian: the boundary blocks ----
% Every second-derivative block is  (1/pi) \int_0^inf t*Im[W(t) phi e^{-i t x0}] dt
% for a block-specific weight W built from M^{-1} and mu_tilde (the extra factor
% it beyond the gradient makes these t-weighted -- hence the s=0 caveat above).
% One array-valued integration returns all blocks; they are then unpacked.
raw_hess=[];
if want_hess
    if s0
        % Robust route: each block is the eigenbasis engine at p0=2. The mode
        % index maps: a1,a2,b1,b2 free indices tie to the
        % mode variables of the M^{-1}/mu_tilde factors in W_{Q2 Q2} etc.
        q0q0=Dterm([],2);                                              % (it)^2 * 1
        q0q1=reshape(evalblock(mk(1,1,true),1,2),d,nx);               % (it)^2 mu_tilde
        % q0q2 == q1q1 == d_q0(d_Q2 F): (it)^2 (M^{-1}+mu_tilde mu_tilde')
        Pblk=[mk(1,[1 1],false),mk(1,[1 2],[true true])];
        q0q2=reshape(evalblock(Pblk,2,2),d,d,nx);
        q1q1=q0q2;
        % q1q2: 2 M^{-1}(a,b)mu_tilde(c) + mu_tilde(a)(M^{-1}+mu_tilde mu_tilde')(b,c)
        q1q2blk=[mk(2,[1 1 2],[false true]), ...
                 mk(1,[1 2 2],[true false]), ...
                 mk(1,[1 2 3],[true true true])];
        q1q2=reshape(evalblock(q1q2blk,3,2),d,d,d,nx);
        % q2q2: the six monomials of W_{Q2 Q2}(a1,a2,b1,b2)
        q2q2blk=[mk(2,[2 1 1 2],[false false]), ...            % 2 Minv(a2,b1)Minv(b2,a1)
                 mk(2,[3 1 1 2],[false true true]), ...        % 2 Minv(a2,b1)mm(b2,a1)
                 mk(2,[3 1 2 3],[true true false]), ...        % 2 mm(a2,b1)Minv(b2,a1)
                 mk(1,[1 1 2 2],[false false]), ...            % Minv(a2,a1)Minv(b2,b1)
                 mk(1,[1 1 3 2],[false true true]), ...        % Minv(a2,a1)mm(b2,b1)
                 mk(1,[2 1 3 3],[true true false]), ...        % mm(a2,a1)Minv(b2,b1)
                 mk(1,[2 1 4 3],[true true true true])];       % mm(a2,a1)mm(b2,b1)
        q2q2=reshape(evalblock(q2q2blk,4,2),d,d,d,d,nx);
    else
        % Inversion route: every block is (1/pi)\int t*Im[W(t) phi e^{-i t x0}] dt
        % for a block weight W built from M^{-1} and mu_tilde (the extra factor
        % it beyond the gradient makes these t-weighted). One array-valued
        % integration returns all blocks; they are then unpacked.
        if ~strcmpi(precision,'basic')
            warning('cdf_grad_norm_quad:hessvpa',...
                'The Hessian uses the ''basic'' integration path; ''precision'' applies to the gradient only.');
        end
        Hraw=integral(@(t) hess_integrand(t,s),0,inf,'ArrayValued',true,'AbsTol',AbsTol,'RelTol',RelTol);
        Hraw=reshape(Hraw,Nt,nx)/pi;
        off=0;
        q0q0=Hraw(1,:);                              off=1;
        q0q1=Hraw(off+(1:d),:);                      off=off+d;
        q0q2=reshape(Hraw(off+(1:d^2),:),d,d,nx);    off=off+d^2;
        q1q1=reshape(Hraw(off+(1:d^2),:),d,d,nx);    off=off+d^2;
        q1q2=reshape(Hraw(off+(1:d^3),:),d,d,d,nx);  off=off+d^3;
        q2q2=reshape(Hraw(off+(1:d^4),:),d,d,d,d,nx);
    end
    q0q2=symm_pages(q0q2);                           % symmetrize the plain matrices
    q1q1=symm_pages(q1q1);
    raw_hess=struct('q0q0',q0q0,'q0q1',q0q1,'q0q2',q0q2,'q1q1',q1q1,'q1q2',q1q2,'q2q2',q2q2);
end

inversion_grad_hess=@inversion_grad_hess_fn;

    % ---- robust s=0 engine --------------------------------------------------
    % A "block" is a list of monomials in M^{-1} and mu_tilde. Each monomial mo
    % records, per free tensor index, which mode variable supplies its u-column
    % (mo.mv), and, per mode variable, whether that factor carries the linear
    % term c=alpha+i t beta (mo.cv true for a mu_tilde factor, false for a bare
    % M^{-1}). evalblock sums over mode assignments; for each it expands the
    % c-factors (cexp) and reads off shifted-dof density derivatives (Dterm).

    function mo=mk(pref,mv,cv)
        mo=struct('pref',pref,'mv',mv,'cv',logical(cv));
    end

    function T=evalblock(monos,F,p0)
        % Assemble a block with F free indices, base it-power p0 (1 gradient,
        % 2 Hessian). Returns [d^F, nx]; caller reshapes to [d,..,d,nx].
        T=zeros(d^F,nx);
        for im=1:numel(monos)
            mo=monos(im); nmv=max(mo.mv);
            for lin=0:(d^nmv-1)                 % odometer over mode assignments
                assign=zeros(1,nmv); r=lin;
                for vv=1:nmv, assign(vv)=mod(r,d)+1; r=floor(r/d); end
                sc=cexp(assign(mo.cv),assign,p0);   % 1-by-nx
                if ~any(sc), continue; end
                rank1=1;                            % rank-1 coeff tensor (u-columns)
                for f=1:F
                    col=U(:,assign(mo.mv(f)));
                    if f==1, rank1=col; else, rank1=kron(col,rank1); end
                end
                T=T+mo.pref*(rank1*sc);
            end
        end
    end

    function val=cexp(cmodes,gmodes,p0)
        % T[(it)^p0 * prod_g(gmodes) * prod_{i in cmodes}(alpha_i+i t beta_i) phi]
        % = sum over subsets S of cmodes: (prod_S beta)(prod_rest alpha)
        %   * Dterm(gmodes, p0+|S|).  gmodes fixes the dof shift; each beta pick
        % raises the it-power (hence the density-derivative order) by one.
        nc=numel(cmodes); val=zeros(1,nx);
        for mask=0:(2^nc-1)
            coef=1; nb=0;
            for ii=1:nc
                if bitget(mask,ii), coef=coef*bet(cmodes(ii)); nb=nb+1;
                else,               coef=coef*alph(cmodes(ii)); end
            end
            if coef~=0, val=val+coef*Dterm(gmodes,p0+nb); end
        end
    end

    function val=Dterm(gmodes,p)
        % T[(it)^p prod_g(gmodes) phi] = (-1)^p f^{(p-1)}_[shift](x), where each
        % g_j advances its component's dof by 2 (rule R1) and (it)^p is p
        % argument-derivatives (rule R2). Zero modes (g_j=1) add no shift.
        bumpvec=zeros(1,numel(k));
        for mm=gmodes
            c=compj(mm); if c>0, bumpvec(c)=bumpvec(c)+2; end
        end
        val=((-1)^p)*fder(bumpvec,p-1);
    end

    function val=fder(bumpvec,n)
        % memoized robust n-th density derivative of the gx2 with k+bumpvec dof
        key=sprintf('%d_',[bumpvec n]);
        if isKey(memo,key), val=memo(key); return; end
        val=reshape(gx2_dens_deriv(x,w,k+bumpvec,l,s,m,n,densopts{:}),1,nx);
        memo(key)=val;
    end

    function A3=symm_pages(A3)
        % symmetrize each d-by-d page (kills round-off in the plain matrix blocks)
        for pp=1:size(A3,3), A3(:,:,pp)=0.5*(A3(:,:,pp)+A3(:,:,pp).'); end
    end

    function out=integrand(t,s_eff)
        % base complex block at scalar t, times phi(t), broadcast over the
        % thresholds x; returns a real [d, d+1, nx] array. Parametrized by
        % s_eff so the same code serves the normal s~=0 case (s_eff=s) and
        % the cusp probe (s_eff -> tiny artificial values, s itself being 0
        % there).
        [Minv,mut,phi]=weights(t,s_eff);
        block=[mut, Minv+mut*mut.'];      % d-by-(d+1): [q1-weight | Q2-weight]
        Bphi=block*phi;
        kern=reshape(exp(-1i*t*x),1,1,nx);
        out=real(Bphi.*kern);
    end

    function [Minv,mut,phi]=weights(t,s_eff)
        % per-node characteristic function and the tilted covariance/mean, from
        % the shared eigen-structure: M^{-1}=S*V*diag(1/(1-2i t d_j))*V'*S.
        phi=gx2char(t,w,k,l,s_eff,m);
        g=1./(1-2i*t*dvals);              % 1-by-d
        Minv=S*(V.*g)*V.'*S;
        pv=SigInv_mu+1i*t*q1c;            % p(t), d-by-1
        mut=Minv*pv;                      % mu_tilde
    end

    function outv=hess_integrand(t,s_eff)
        [Minv,mut,phi]=weights(t,s_eff);
        P=Minv+mut*mut.';                 % M^{-1}+mu_tilde*mu_tilde'
        Pmm=mut*mut.';
        Wq1q2=zeros(d,d,d);
        for a=1:d, for b=1:d, for c=1:d
            Wq1q2(a,b,c)=2*Minv(a,b)*mut(c)+mut(a)*P(b,c);
        end, end, end
        WQ2=zeros(d,d,d,d);
        for a1=1:d, for a2=1:d, for b1=1:d, for b2=1:d
            WQ2(a1,a2,b1,b2)=2*(Minv(a2,b1)*Minv(b2,a1)+Minv(a2,b1)*Pmm(b2,a1) ...
                +Pmm(a2,b1)*Minv(b2,a1))+P(a2,a1)*P(b2,b1);
        end, end, end, end
        Wvec=[1; mut(:); P(:); P(:); Wq1q2(:); WQ2(:)];   % Nt-by-1 complex
        outv=zeros(numel(Wvec),nx);
        for ix=1:nx
            outv(:,ix)=t*imag(Wvec*phi*exp(-1i*t*x(ix)));
        end
    end

    function [gr,hr]=inversion_grad_hess_fn(s_eff,want_hess_probe)
        % Full grad (and, if want_hess_probe, Hessian) via the s~=0 inversion
        % route at an arbitrary s_eff, in the same raw (pre-squeeze,
        % trailing-nx) shape returned above. Used both for the normal s~=0
        % case (s_eff=s, called inline above rather than through here) and,
        % when at_cusp, as the cusp-convergence probe (via resolve_cusp_terms).
        g0r=reshape(-gx2pdf(x,w,k,l,s_eff,m,opts{:}),1,nx);
        A=integral(@(t) integrand(t,s_eff),0,inf,'ArrayValued',true,'AbsTol',AbsTol,'RelTol',RelTol);
        A=-A/pi;
        Gq1r=reshape(A(:,1,:),d,nx);
        Gq2r=symm_pages(reshape(A(:,2:d+1,:),d,d,nx));
        gr=struct('q0',g0r,'q1',Gq1r,'q2',Gq2r);
        if ~want_hess_probe, hr=[]; return; end
        Hraw=integral(@(t) hess_integrand(t,s_eff),0,inf,'ArrayValued',true,'AbsTol',AbsTol,'RelTol',RelTol);
        Hraw=reshape(Hraw,Nt,nx)/pi;
        off=0;
        q0q0r=Hraw(1,:);                              off=1;
        q0q1r=Hraw(off+(1:d),:);                      off=off+d;
        q0q2r=reshape(Hraw(off+(1:d^2),:),d,d,nx);    off=off+d^2;
        q1q1r=reshape(Hraw(off+(1:d^2),:),d,d,nx);    off=off+d^2;
        q1q2r=reshape(Hraw(off+(1:d^3),:),d,d,d,nx);  off=off+d^3;
        q2q2r=reshape(Hraw(off+(1:d^4),:),d,d,d,d,nx);
        q0q2r=symm_pages(q0q2r);
        q1q1r=symm_pages(q1q1r);
        hr=struct('q0q0',q0q0r,'q0q1',q0q1r,'q0q2',q0q2r,'q1q1',q1q1r,'q1q2',q1q2r,'q2q2',q2q2r);
    end

    function A=vpa_integrate()
        % variable-precision path: vpaintegral has no array-valued mode, so
        % integrate each entry (and threshold) separately. The function-handle
        % form infers the integration variable, so we avoid a `syms` here
        % (which a nested function's static workspace would reject).
        A=zeros(d,d+1,nx);
        for ix=1:nx
            xi=x(ix);
            for r=1:d
                for c=1:d+1
                    A(r,c,ix)=double(vpaintegral(@(t) entry_real(t,r,c,xi),...
                        0,inf,'AbsTol',AbsTol,'RelTol',RelTol));
                end
            end
        end
    end

    function y=entry_real(t,r,c,xi)
        Bphi=block_phi(t);
        y=real(Bphi(r,c).*exp(-1i*t*xi));
    end
end
