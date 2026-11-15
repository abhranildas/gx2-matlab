function f=gx2_imhof_integrand(u,x,w,k,l,s,m,output,idx,nx)
% Imhof integrand for the generalized chi-square inversion (w, k, l are
% column vectors here). Beyond the plain cdf/pdf it also returns the exact
% integrands for parameter-gradient and -Hessian components (no finite
% differencing).
%
% output:
%   'cdf'      cdf integrand.
%   'pdf'      pdf integrand.
%   'dens'     nx-th x-derivative of the pdf, f^(nx) (nx>=1 gives f',f'',...).
%   'k_deriv'  d/dk_(idx) of the cdf, with nx extra x-derivatives (idx scalar).
%   'kk_deriv' d^2/(dk_(idx1) dk_(idx2)) of the cdf, with nx extra x-derivatives.
% idx: component index (k_deriv) or [i j] (kk_deriv). nx: x-derivative order, default 0.
%
% The derivative modes use the complex integrand Z=exp(i*theta)/rho together
% with the rule (R2) that each x-derivative multiplies Z by -(i*u/2), and that
% d/dk_j brings down the factor ell_j=-1/2*log(1-i*w_j*u).
%
% Every mode is real(H(u).*exp(-1i*c*u)) with c=(x-m)/2 and the
% x-independent amplitude H from gx2_imhof_amplitude.

if nargin<9, idx=[]; end
if nargin<10 || isempty(nx), nx=0; end

H=gx2_imhof_amplitude(u,w,k,l,s,output,idx,nx);
f=real(H.*exp(-1i*u.*(x-m)/2));
end
