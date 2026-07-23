function f=gx2_imhof_integrand(u,x,w,k,lambda,s,m,output,idx,nx)
% Imhof integrand for the generalized chi-square inversion (w, k, lambda are
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

if nargin<9, idx=[]; end
if nargin<10 || isempty(nx), nx=0; end

theta=sum(k.*atan(w*u)+(lambda.*(w*u))./(1+w.^2*u.^2),1)/2+u*(m-x)/2;
rho=prod(((1+w.^2*u.^2).^(k/4)).*exp(((w.^2*u.^2).*lambda)./(2*(1+w.^2*u.^2))),1) .* exp(u.^2*s^2/8);

if strcmpi(output,'cdf')
    f=sin(theta)./(u.*rho);
elseif strcmpi(output,'pdf')
    f=cos(theta)./rho;
else
    Z=exp(1i*theta)./rho;          % phi(t)*exp(-i*t*x), with u=2t
    dx=(-(1i*u/2)).^nx;            % nx x-derivatives (R2)
    if strcmpi(output,'dens')
        f=real(dx.*Z);
    elseif strcmpi(output,'k_deriv')
        ell=-0.5*log(1-1i*w(idx)*u);
        f=-imag(ell.*dx.*Z)./u;
    elseif strcmpi(output,'kk_deriv')
        ell=-0.5*log(1-1i*w(idx(1))*u).*(-0.5*log(1-1i*w(idx(2))*u));
        f=-imag(ell.*dx.*Z)./u;
    end
end
end
