function fd=gx2_dens_deriv(x,w,k,l,s,m,nx,varargin)

% GX2_DENS_DERIV Robust nx-th derivative in x of the generalized chi-square pdf.
%
% Abhranil Das
% Center for Perceptual Systems, University of Texas at Austin
% Comments, questions, bugs to abhranil.das@utexas.edu
%
% Single entry point for the density x-derivatives f', f'', f''' used by the
% Hessian routines. It picks between two exact methods:
%
%   - the Gil-Pelaez (Imhof) t-weighted inversion, used whenever it converges
%     comfortably: s~=0 (Gaussian damping), or s=0 with total dof large
%     relative to the derivative order;
%   - the differentiated shifted-dof (Ruben) series, used at s=0 with small
%     total dof, where the inversion integrand loses convergence. For
%     mixed-sign weights the variable is split as q = q_+ - q_-, each part a
%     same-sign (elliptical) gx2 that Ruben's series handles, and the density
%     is their cross-correlation with the derivatives falling on q_+.
%
% Usage:
% fd=gx2_dens_deriv(x,w,k,l,s,m,nx)
% fd=gx2_dens_deriv(x,w,k,l,s,m,nx,'AbsTol',0,'RelTol',1e-9,'precision','vpa')
%
% Inputs mirror gx2pdf, plus nx (the derivative order, 0 gives the pdf).
% Optional 'AbsTol','RelTol','precision' pass through to the Imhof route and
% the mixed-sign convolution quadrature.

parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'x',@(x) isreal(x));
addRequired(parser,'w',@(x) isreal(x) && isrow(x));
addRequired(parser,'k',@(x) isreal(x) && isrow(x));
addRequired(parser,'l',@(x) isreal(x) && isrow(x));
addRequired(parser,'s',@(x) isreal(x) && isscalar(x));
addRequired(parser,'m',@(x) isreal(x) && isscalar(x));
addRequired(parser,'nx',@(x) isscalar(x) && (x>=0) && (x==round(x)));
addParameter(parser,'AbsTol',1e-10,@(x) isreal(x) && isscalar(x) && (x>=0));
addParameter(parser,'RelTol',1e-6,@(x) isreal(x) && isscalar(x) && (x>=0));
addParameter(parser,'n_ruben',1e3,@(x) isscalar(x) && (x>0) && (x==round(x)));
parse(parser,x,w,k,l,s,m,nx,varargin{:});
AbsTol=parser.Results.AbsTol;
RelTol=parser.Results.RelTol;
n_ruben=parser.Results.n_ruben;

D=sum(k);   % total degrees of freedom

% Method choice at s=0 (a nonzero s adds Gaussian damping that makes the
% inversion converge fast for every order, so the inversion is used then).
% The Imhof density-derivative integrand behaves like u^{nx-D/2} as u->inf, so
% at s=0 it is only conditionally convergent for D<=2*nx+2 and slowly
% convergent just above; the differentiated series is needed there. For
% same-sign (elliptical) weights the series applies at ANY dof and is 1-3
% orders of magnitude faster than the inversion (which converges but crawls
% for a slowly-decaying integrand), so prefer it regardless of dof. For mixed
% signs the series needs a convolution, so fall back to it only in the
% small-dof regime where the inversion loses convergence.
same_sign = all(w>0) || all(w<0);
use_series = (s==0) && (same_sign || (D <= 2*nx+3));

if ~use_series
    if nx==0
        % the plain pdf: use gx2pdf's default dispatch (exact where possible)
        fd=gx2pdf(x,w,k,l,s,m,varargin{:});
    else
        fd=gx2_imhof(x,w,k,l,s,m,varargin{:},'output','dens','nx',nx);
    end
    return;
end

% ---- s=0 series route ----
pos=w>0; neg=w<0;
if all(pos) || all(neg)
    % same-sign (elliptical): differentiate Ruben's series directly
    fd=gx2_ruben(x,w,k,l,m,'output','pdf','nx',nx,'n_ruben',n_ruben);
else
    % mixed sign: q = q_+ - q_-, where q_+ = (positive-weight part) + m has
    % support [m,inf) and q_- = (negated negative-weight part) has support
    % [0,inf). The density is their cross-correlation, whose x-derivatives may
    % be carried on either part. A low-dof density has a non-integrable edge
    % singularity once differentiated, so we carry the derivatives on whichever
    % part is sampled in its smooth interior, away from its own edge (the
    % differentiated factor then never meets that singularity). For a threshold
    % x below the q_+ floor m we differentiate q_-, otherwise q_+:
    %   x <  m:  f^(nx)(x) = (-1)^nx int f_{q+}(x+v) f_{q-}^(nx)(v) dv
    %   x >= m:  f^(nx)(x) =          int f_{q+}^(nx)(x+v) f_{q-}(v) dv
    % (see section 1.5). At x=m exactly, both floors align into a genuine cusp;
    % it is measure-zero and not hit in practice.
    wp=w(pos);  kp=k(pos);  lp=l(pos);
    wn=-w(neg); kn=k(neg);  ln=l(neg);   % negate -> positive weights

    % Ruben's series coefficients depend only on (w,k,l), not on the
    % evaluation point -- compute them once per gx2_dens_deriv call and reuse
    % across every scalar quadrature callback below, rather than rebuilding
    % the series from scratch on each of the (potentially hundreds of) points
    % the inner integral() sweeps per integral.
    coeffs_p=gx2_ruben_coeffs(wp,kp,lp,'n_ruben',n_ruben);
    coeffs_n=gx2_ruben_coeffs(wn,kn,ln,'n_ruben',n_ruben);
    fp=@(y,n) gx2_ruben_eval(coeffs_p,y,m,'output','pdf','nx',n);
    fm=@(v,n) gx2_ruben_eval(coeffs_n,v,0,'output','pdf','nx',n);
    fd=arrayfun(@(xx) conv_dens_deriv(xx,m,nx,fp,fm,AbsTol,RelTol),x);
end
fd=reshape(fd,size(x));
end

function val=conv_dens_deriv(xx,m,nx,fp,fm,AbsTol,RelTol)
% One point of the mixed-sign cross-correlation, with the nx x-derivatives
% carried on the part sampled in its interior (see the parent function).
if xx<m
    % Differentiate q_-. Since f_{q+}(xx+v)=0 for v<m-xx, start the integral at
    % that floor; the (integrable) q_+ edge is then the lower endpoint, and the
    % singular f_{q-}^(nx) is evaluated only for v>=m-xx>0, in q_-'s interior.
    val=(-1)^nx*integral(@(v) fp(xx+v,0).*fm(v,nx),m-xx,inf,...
        'AbsTol',AbsTol,'RelTol',RelTol);
else
    % Differentiate q_+, which is sampled strictly inside its support (xx+v>=m).
    val=integral(@(v) fp(xx+v,nx).*fm(v,0),0,inf,...
        'AbsTol',AbsTol,'RelTol',RelTol);
end
end
