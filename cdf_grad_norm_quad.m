function [grad,hess]=cdf_grad_norm_quad(x,mu,v,quad,varargin)

% CDF_GRAD_NORM_QUAD Gradient (and optionally Hessian) of the cdf of a quadratic
% form q(x)=x'*Q2*x + q1'*x + q0 of a normal vector x~N(mu,v), with respect to
% the quadratic's coefficients Q2, q1, q0 (holding mu and v fixed).
%
% F(x0)=P(q(x)<=x0) is the probability content of the normal in the quadratic
% region q(x)<=x0. This returns its derivatives with respect to Q2, q1 and q0,
% computed exactly (no finite differencing).
%
% Abhranil Das
% Center for Perceptual Systems, University of Texas at Austin
% Comments, questions, bugs to abhranil.das@utexas.edu
% If you use this code, please cite:
% 1. <a href="matlab:web('https://doi.org/10.1167/jov.21.10.1','-browser')"
% >A method to integrate and classify normal distributions</a>
% 2. <a href="matlab:web('https://www.tandfonline.com/doi/abs/10.1080/00949655.2025.2501401','-browser')"
% >New methods to compute the generalized chi-square distribution</a>
%
% Usage:
% grad=cdf_grad_norm_quad(x,mu,v,quad)
% grad=cdf_grad_norm_quad(x,mu,v,quad,'wrt',{'q2','q0'})
% grad=cdf_grad_norm_quad(x,mu,v,quad,'AbsTol',0,'RelTol',1e-7,'precision','vpa')
%
% Example:
% mu=[1;2]; v=[2 1; 1 3];
% quad.q2=[1 1; 1 1]; quad.q1=[-1;0]; quad.q0=-1;
% grad=cdf_grad_norm_quad(0,mu,v,quad)
%
% Required inputs:
% x         array of thresholds x0 at which to evaluate the gradient of the cdf
% mu        column vector of normal mean
% v         normal covariance matrix
% quad      struct with quadratic form coefficients:
%               q2      matrix of quadratic coefficients (symmetrized internally)
%               q1      column vector of linear coefficients
%               q0      scalar constant
%
% Optional name-value inputs:
% wrt       cell array selecting which coefficient groups to differentiate
%           with respect to, drawn from {'q2','q1','q0'}. Default is all.
% AbsTol    absolute error tolerance for the underlying integrals. Default=1e-10.
% RelTol    relative error tolerance for the underlying integrals. Default=1e-6.
%           The absolute OR the relative tolerance is satisfied.
% precision 'basic' (default) uses double precision, 'vpa' uses variable precision.
%
% Outputs:
% grad      struct mirroring quad, holding the cdf gradient:
%               q2      symmetric d-by-d matrix G = dF/dQ2, in the sense
%                       dF ~ trace(G*dQ2) for symmetric perturbations dQ2
%                       (so a lone off-diagonal (Q2)_ab sees 2*G_ab).
%               q1      d-vector dF/dq1
%               q0      scalar dF/dq0 = -pdf of q(x) at the threshold
%           For an array x each field carries a trailing dimension over the
%           points (q2: [d,d,numel(x)]; q1: [d,numel(x)]; q0: [1,numel(x)]).
%           Groups omitted by 'wrt' are absent from the struct.
%
% See also:
% cdf_grad_gx2, norm_quad_to_gx2_params, gx2cdf, gx2pdf, gx2char

parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'x',@(x) isreal(x));
addRequired(parser,'mu',@isnumeric);
addRequired(parser,'v',@isnumeric);
addRequired(parser,'quad',@isstruct);
groups={'q2','q1','q0'};
addParameter(parser,'wrt',groups,@(c) iscell(c) && all(ismember(lower(c),groups)));
addParameter(parser,'AbsTol',1e-10,@(x) isreal(x) && isscalar(x) && (x>=0));
addParameter(parser,'RelTol',1e-6,@(x) isreal(x) && isscalar(x) && (x>=0));
addParameter(parser,'precision','basic',@(x) strcmpi(x,'basic')||strcmpi(x,'vpa'));
parse(parser,x,mu,v,quad,varargin{:});

wrt=parser.Results.wrt;
AbsTol=parser.Results.AbsTol;
RelTol=parser.Results.RelTol;
precision=parser.Results.precision;

mu=mu(:);
q1c=quad.q1(:);
x=x(:)';                 % row of thresholds
nx=numel(x);
d=numel(mu);
wanted=@(g) any(strcmpi(g,wrt));
opts={'AbsTol',AbsTol,'RelTol',RelTol,'precision',precision};

% convert to gx2 params, and reuse the standardized-quadratic eigen-structure
% (S=Sigma^{1/2}, V, and the full eigenvalues d of S*Q2*S) for the per-node
% M^{-1}(t)=S*V*diag(1/(1-2i*t*d_j))*V'*S -- no d-by-d inverse per node.
[w,k,lambda,s,m,aux]=norm_quad_to_gx2_params(mu,v,quad);
S=aux.S; V=aux.V; dvals=aux.d(:).';     % dvals is 1-by-d
SigInv_mu=v\mu;                          % Sigma^{-1}*mu (real d-vector)

grad=struct();

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
    densopts={'AbsTol',AbsTol,'RelTol',RelTol};
end

% q0 block: dF/dq0 = -f(x0). Since q0 shifts q rigidly, this is just -pdf.
if wanted('q0')
    if s0, g0=Dterm([],1); else, g0=reshape(-gx2pdf(x,w,k,lambda,s,m,opts{:}),1,nx); end
    grad.q0=g0;
    if nx==1, grad.q0=g0(1); end
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
            Gq1=reshape(evalblock(mk(1,1,true),1,1),d,nx);
            if nx==1, grad.q1=Gq1(:,1); else, grad.q1=Gq1; end
        end
        if wanted('q2')
            Gq2=symm_pages(reshape(evalblock([mk(1,[1 1],false),mk(1,[1 2],[true true])],2,1),d,d,nx));
            if nx==1, grad.q2=Gq2(:,:,1); else, grad.q2=Gq2; end
        end
    else
        if strcmpi(precision,'basic')
            A=integral(@integrand,0,inf,'ArrayValued',true,'AbsTol',AbsTol,'RelTol',RelTol);
        else
            A=vpa_integrate();
        end
        A=-A/pi;                              % [d, d+1, nx]
        if wanted('q1')
            Gq1=reshape(A(:,1,:),d,nx);
            if nx==1, grad.q1=Gq1(:,1); else, grad.q1=Gq1; end
        end
        if wanted('q2')
            Gq2=symm_pages(reshape(A(:,2:d+1,:),d,d,nx));
            if nx==1, grad.q2=Gq2(:,:,1); else, grad.q2=Gq2; end
        end
    end
end

% ---- Hessian (2nd output): the boundary blocks ----
% Every second-derivative block is  (1/pi) \int_0^inf t*Im[W(t) phi e^{-i t x0}] dt
% for a block-specific weight W built from M^{-1} and mu_tilde (the extra factor
% it beyond the gradient makes these t-weighted -- hence the s=0 caveat above).
% One array-valued integration returns all blocks; they are then unpacked.
if nargout>=2
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
        Nt=1+d+d^2+d^2+d^3+d^4;          % q0q0, q0q1, q0q2, q1q1, q1q2, q2q2
        if ~strcmpi(precision,'basic')
            warning('cdf_grad_norm_quad:hessvpa',...
                'The Hessian uses the ''basic'' integration path; ''precision'' applies to the gradient only.');
        end
        Hraw=integral(@hess_integrand,0,inf,'ArrayValued',true,'AbsTol',AbsTol,'RelTol',RelTol);
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
    hess=struct('q0q0',q0q0,'q0q1',q0q1,'q0q2',q0q2,'q1q1',q1q1,'q1q2',q1q2,'q2q2',q2q2);
    if nx==1                                     % drop the trailing singleton
        hess.q0q0=q0q0(1); hess.q0q1=q0q1(:,1);
        hess.q0q2=q0q2(:,:,1); hess.q1q1=q1q1(:,:,1);
        hess.q1q2=q1q2(:,:,:,1); hess.q2q2=q2q2(:,:,:,:,1);
    end
end

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
        val=reshape(gx2_dens_deriv(x,w,k+bumpvec,lambda,s,m,n,densopts{:}),1,nx);
        memo(key)=val;
    end

    function A3=symm_pages(A3)
        % symmetrize each d-by-d page (kills round-off in the plain matrix blocks)
        for pp=1:size(A3,3), A3(:,:,pp)=0.5*(A3(:,:,pp)+A3(:,:,pp).'); end
    end

    function out=integrand(t)
        % base complex block at scalar t, times phi(t), broadcast over the
        % thresholds x; returns a real [d, d+1, nx] array.
        [Minv,mut,phi]=weights(t);
        block=[mut, Minv+mut*mut.'];      % d-by-(d+1): [q1-weight | Q2-weight]
        Bphi=block*phi;
        kern=reshape(exp(-1i*t*x),1,1,nx);
        out=real(Bphi.*kern);
    end

    function [Minv,mut,phi]=weights(t)
        % per-node characteristic function and the tilted covariance/mean, from
        % the shared eigen-structure: M^{-1}=S*V*diag(1/(1-2i t d_j))*V'*S.
        phi=gx2char(t,w,k,lambda,s,m);
        g=1./(1-2i*t*dvals);              % 1-by-d
        Minv=S*(V.*g)*V.'*S;
        pv=SigInv_mu+1i*t*q1c;            % p(t), d-by-1
        mut=Minv*pv;                      % mu_tilde
    end

    function outv=hess_integrand(t)
        [Minv,mut,phi]=weights(t);
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
