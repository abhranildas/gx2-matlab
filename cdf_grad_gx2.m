function [grad,hess]=cdf_grad_gx2(x,w,k,lambda,s,m,varargin)

% CDF_GRAD_GX2 Returns the gradient (and optionally the Hessian) of the cdf of
% a generalized chi-squared distribution with respect to its parameters w, k,
% lambda, s and m. These are computed exactly, with no finite differencing.
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
% grad=cdf_grad_gx2(x,w,k,lambda,s,m)
% [grad,hess]=cdf_grad_gx2(x,w,k,lambda,s,m)
% grad=cdf_grad_gx2(x,w,k,lambda,s,m,'wrt',{'s','m'})
% grad=cdf_grad_gx2(x,w,k,lambda,s,m,'AbsTol',0,'RelTol',1e-7,'precision','vpa')
%
% Example:
% [grad,hess]=cdf_grad_gx2(25,[1 -5 2],[1 2 3],[2 3 7],0,5)
%
% Required inputs:
% x         array of points at which to evaluate the gradient/Hessian of the cdf
% w         row vector of weights of the non-central chi-squares
% k         row vector of degrees of freedom of the non-central chi-squares
% lambda    row vector of non-centrality paramaters (sum of squares of
%           means) of the non-central chi-squares
% s         scale of normal term
% m         offset
%
% Optional name-value inputs:
% wrt       cell array selecting which parameter groups to differentiate
%           with respect to, drawn from {'w','k','lambda','s','m'}.
%           Default is all of them. Only the requested groups are returned,
%           in the canonical order below with the unrequested groups omitted;
%           the Hessian is the corresponding principal submatrix. Use this to
%           avoid returning the (possibly many) per-component derivatives when
%           only a few are wanted.
% AbsTol    absolute error tolerance for the underlying integrals. Default=1e-10.
% RelTol    relative error tolerance for the underlying integrals. Default=1e-6.
%           The absolute OR the relative tolerance is satisfied.
% precision 'basic' (default) uses double precision, 'vpa' uses variable precision.
%
% Outputs:
% grad      gradient of the cdf, as a flat numeric array. For a scalar x it
%           is a column vector; for an array x it has one column per point,
%           i.e. size [P, numel(x)], where P is the number of requested
%           parameters. The rows are stacked in the canonical order
%
%               [ dF/dw_1 ... dF/dw_n,
%                 dF/dk_1 ... dF/dk_n,
%                 dF/dlambda_1 ... dF/dlambda_n,
%                 dF/ds,
%                 dF/dm ]
%
%           where n=numel(w), for a total length 3n+2 when all groups are
%           requested. When 'wrt' omits a group, its rows are dropped and
%           the remaining rows keep this relative order.
% hess      (optional second output) Hessian of the cdf: the symmetric matrix
%           of second derivatives d^2F/(da db) over the same parameters and in
%           the same canonical order as grad. For a scalar x it is [P, P]; for
%           an array x it is [P, P, numel(x)]. Requesting hess computes all
%           blocks regardless of 'wrt', then returns the requested submatrix.
%
% See also:
% gx2cdf, gx2pdf, gx2char

parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'x',@(x) isreal(x));
addRequired(parser,'w',@(x) isreal(x) && isrow(x));
addRequired(parser,'k',@(x) isreal(x) && isrow(x));
addRequired(parser,'lambda',@(x) isreal(x) && isrow(x));
addRequired(parser,'s',@(x) isreal(x) && isscalar(x));
addRequired(parser,'m',@(x) isreal(x) && isscalar(x));
groups={'w','k','lambda','s','m'};
addParameter(parser,'wrt',groups,@(c) iscell(c) && all(ismember(lower(c),groups)));
addParameter(parser,'AbsTol',1e-10,@(x) isreal(x) && isscalar(x) && (x>=0));
addParameter(parser,'RelTol',1e-6,@(x) isreal(x) && isscalar(x) && (x>=0));
addParameter(parser,'precision','basic',@(x) strcmpi(x,'basic')||strcmpi(x,'vpa'));
parse(parser,x,w,k,lambda,s,m,varargin{:});

wrt=parser.Results.wrt;
AbsTol=parser.Results.AbsTol;
RelTol=parser.Results.RelTol;
precision=parser.Results.precision;

x=x(:)';        % work with a row of evaluation points
nx=numel(x);
n=numel(w);
wanted=@(g) any(strcmpi(g,wrt));

% Options forwarded to the underlying calls. The cdf/pdf building blocks
% (F, f, and their shifted-dof versions) are plain gx2 cdf/pdf values, so
% they use gx2cdf/gx2pdf's default 'auto' dispatch, which picks exact closed
% forms (ncx2, normal) or Ruben's series where available and Imhof otherwise
% -- all deterministic, so no Monte-Carlo noise enters the finite differences
% in dof. Only the special density/k-derivative integrands are intrinsically
% Imhof, and call gx2_imhof directly.
opts={'AbsTol',AbsTol,'RelTol',RelTol,'precision',precision};
imhof_opts=opts;

% ---------------------------- gradient (1st output) ----------------------------
% base blocks, computed only if a group needs them
F=[]; f=[];
if wanted('lambda'), F=gx2cdf(x,w,k,lambda,s,m,opts{:}); end     % dF/dlambda needs F
if wanted('m'),      f=gx2_dens_deriv(x,w,k,lambda,s,m,0,opts{:}); end   % dF/dm = -f (robust pdf)

% dF/dw_j = -k_j f_[k_j+2] - lambda_j f_[k_j+4]
if wanted('w')
    gw=zeros(n,nx);
    for j=1:n
        kp2=k; kp2(j)=kp2(j)+2;
        kp4=k; kp4(j)=kp4(j)+4;
        f2=gx2_dens_deriv(x,w,kp2,lambda,s,m,0,opts{:});
        f4=gx2_dens_deriv(x,w,kp4,lambda,s,m,0,opts{:});
        gw(j,:)=-k(j)*f2-lambda(j)*f4;
    end
end

% dF/dk_j = convergent Imhof integral (no closed shift form)
if wanted('k')
    gk=zeros(n,nx);
    for j=1:n
        gk(j,:)=gx2_imhof(x,w,k,lambda,s,m,imhof_opts{:},'output','k_deriv','idx',j);
    end
end

% dF/dlambda_j = (F_[k_j+2] - F)/2
if wanted('lambda')
    gl=zeros(n,nx);
    for j=1:n
        kp2=k; kp2(j)=kp2(j)+2;
        F2=gx2cdf(x,w,kp2,lambda,s,m,opts{:});
        gl(j,:)=0.5*(F2-F);
    end
end

% dF/ds = s f', with f' the density x-derivative from a t-weighted integrand.
% dF/ds is identically 0 at s=0, so skip f' there (also its worst-conditioned
% point).
if wanted('s')
    if s==0
        gs=zeros(1,nx);
    else
        fprime=gx2_dens_deriv(x,w,k,lambda,s,m,1,imhof_opts{:});
        gs=s*fprime;
    end
end

% dF/dm = -f
if wanted('m')
    gm=-f;
end

% stack in canonical order, dropping omitted groups; build the matching index
% list 'sel' into the full [w;k;lambda;s;m] ordering for the Hessian subset
grad=[]; sel=[];
if wanted('w'),      grad=[grad; gw]; sel=[sel, 1:n];        end
if wanted('k'),      grad=[grad; gk]; sel=[sel, n+(1:n)];    end
if wanted('lambda'), grad=[grad; gl]; sel=[sel, 2*n+(1:n)];  end
if wanted('s'),      grad=[grad; gs]; sel=[sel, 3*n+1];      end
if wanted('m'),      grad=[grad; gm]; sel=[sel, 3*n+2];      end

% ---------------------------- Hessian (2nd output) -----------------------------
if nargout>=2
    P0=3*n+2;
    H=zeros(P0,P0,nx);
    IW=@(j)j; IK=@(j)n+j; IL=@(j)2*n+j; IS=3*n+1; IM=3*n+2;   % full-order indices
    sh=@(kk,j,d) kk+d*((1:n)==j);                            % add d to k(j)

    % building blocks (each returns a 1 x nx row over the evaluation points)
    Fh   =@(kk)    gx2cdf(x,w,kk,lambda,s,m,opts{:});
    fh   =@(kk)    gx2_dens_deriv(x,w,kk,lambda,s,m,0,opts{:});   % robust pdf
    fp   =@(kk)    gx2_dens_deriv(x,w,kk,lambda,s,m,1,imhof_opts{:});   % robust f'
    fpp  =@(kk)    gx2_dens_deriv(x,w,kk,lambda,s,m,2,imhof_opts{:});   % robust f''
    fppp =@(kk)    gx2_dens_deriv(x,w,kk,lambda,s,m,3,imhof_opts{:});   % robust f'''
    dkF  =@(kk,j)  gx2_imhof(x,w,kk,lambda,s,m,imhof_opts{:},'output','k_deriv','idx',j,'nx',0);
    dkf  =@(kk,j)  gx2_imhof(x,w,kk,lambda,s,m,imhof_opts{:},'output','k_deriv','idx',j,'nx',1);
    dkFxx=@(kk,j)  gx2_imhof(x,w,kk,lambda,s,m,imhof_opts{:},'output','k_deriv','idx',j,'nx',2);
    dkkF =@(kk,i,j)gx2_imhof(x,w,kk,lambda,s,m,imhof_opts{:},'output','kk_deriv','idx',[i j],'nx',0);

    F0=Fh(k); f0=fh(k); fp0=fp(k);

    % global. At s=0 the s-coupled entries vanish (they all carry a factor s),
    % so set them to 0 without evaluating the higher density derivatives f'',
    % f''' -- which are both unneeded there and worst-conditioned at s=0.
    z=zeros(1,nx);
    put(IM,IM, fp0);                                  % H_mm = f'
    if s==0
        put(IS,IM, z);                                % H_ms = -s f'' = 0
        put(IS,IS, fp0);                              % H_ss = f' + s^2 f''' = f'
    else
        put(IS,IM, -s*fpp(k));                        % H_ms = -s f''
        put(IS,IS, fp0+s^2*fppp(k));                  % H_ss = f' + s^2 f'''
    end

    for j=1:n
        kj=k(j); lj=lambda(j);
        kp2=sh(k,j,2); kp4=sh(k,j,4); kp6=sh(k,j,6); kp8=sh(k,j,8);
        % global x component
        put(IM,IL(j), -0.5*(fh(kp2)-f0));                                    % H_m,lambda_j
        put(IM,IW(j), kj*fp(kp2)+lj*fp(kp4));                                % H_m,w_j
        % same component
        put(IL(j),IL(j), 0.25*(Fh(kp4)-2*Fh(kp2)+F0));                       % H_lambda_j,lambda_j
        put(IL(j),IW(j), 0.5*kj*fh(kp2)+0.5*(lj-kj-2)*fh(kp4)-0.5*lj*fh(kp6));% H_lambda_j,w_j
        put(IW(j),IW(j), kj*(kj+2)*fp(kp4)+2*lj*(kj+2)*fp(kp6)+lj^2*fp(kp8)); % H_w_j,w_j
        % k-blocks (global/same component)
        put(IM,IK(j), -dkf(k,j));                                            % H_m,k_j = -d_x d_k F
        put(IL(j),IK(j), 0.5*(dkF(kp2,j)-dkF(k,j)));                         % H_lambda_j,k_j
        put(IW(j),IK(j), -fh(kp2)-kj*dkf(kp2,j)-lj*dkf(kp4,j));              % H_w_j,k_j
        put(IK(j),IK(j), dkkF(k,j,j));                                       % H_k_j,k_j (log^2 weight)
        % s-coupled entries (all carry a factor s -> 0 at s=0)
        if s==0
            put(IS,IL(j), z); put(IS,IW(j), z); put(IS,IK(j), z);
        else
            put(IS,IL(j), 0.5*s*(fp(kp2)-fp0));                              % H_s,lambda_j
            put(IS,IW(j), -s*(kj*fpp(kp2)+lj*fpp(kp4)));                     % H_s,w_j
            put(IS,IK(j), s*dkFxx(k,j));                                     % H_s,k_j = s d_x^2 d_k F
        end
    end

    % cross component (i ~= j)
    for i=1:n
        for j=1:n
            if i==j, continue; end
            ki=k(i); li=lambda(i); kj=k(j); lj=lambda(j);
            kip2=sh(k,i,2); kip4=sh(k,i,4); kjp2=sh(k,j,2); kjp4=sh(k,j,4);
            % non-symmetric-in-(i,j) blocks: fill for every ordered pair
            put(IL(i),IW(j), -0.5*(kj*(fh(sh(kip2,j,2))-fh(kjp2))+lj*(fh(sh(kip2,j,4))-fh(kjp4)))); % H_lambda_i,w_j
            put(IL(i),IK(j), 0.5*(dkF(kip2,j)-dkF(k,j)));                                            % H_lambda_i,k_j
            put(IW(i),IK(j), -ki*dkf(kip2,j)-li*dkf(kip4,j));                                        % H_w_i,k_j
            % symmetric-in-(i,j) blocks: fill once
            if i<j
                put(IL(i),IL(j), 0.25*(Fh(sh(kip2,j,2))-Fh(kip2)-Fh(kjp2)+F0));                      % H_lambda_i,lambda_j
                put(IW(i),IW(j), ki*kj*fp(sh(kip2,j,2))+ki*lj*fp(sh(kip2,j,4)) ...
                               +li*kj*fp(sh(kip4,j,2))+li*lj*fp(sh(kip4,j,4)));                      % H_w_i,w_j
                put(IK(i),IK(j), dkkF(k,i,j));                                                       % H_k_i,k_j
            end
        end
    end

    hess=H(sel,sel,:);
    if nx==1, hess=hess(:,:,1); end
end

    function put(a,b,val)
        % place a symmetric Hessian entry (val is 1 x nx over eval points)
        H(a,b,:)=reshape(val,1,1,[]);
        if a~=b
            H(b,a,:)=reshape(val,1,1,[]);
        end
    end
end
