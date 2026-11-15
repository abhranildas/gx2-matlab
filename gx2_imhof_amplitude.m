function H=gx2_imhof_amplitude(u,w,k,l,s,output,idx,nx)
% Complex amplitude H(u) of the Imhof integrand, which is
% real(H(u).*exp(-1i*c*u)) with c=(x-m)/2. H does not depend on x or m and
% varies slowly in u; all the oscillation of the integrand is in the factor
% exp(-1i*c*u). gx2_imhof uses this split to integrate the oscillatory tail
% cycle by cycle.
%
% w, k, l are column vectors; u is a scalar or row vector (or a symbolic
% scalar, for the vpa path). output, idx, nx are as in gx2_imhof_integrand.
%
% Writing Z=exp(1i*theta0)/rho, with theta0 the x-free part of the phase,
% H=q*Z, where q=-1i/u for the cdf (so real(H.*exp(-1i*c*u)) is
% sin(theta)/(u*rho)), 1 for the pdf, the x-derivative factor (-1i*u/2)^nx
% for 'dens', and 1i*ell*(-1i*u/2)^nx/u for the k-derivatives
% (real(1i*y)=-imag(y)).
%
% rho is computed in log form, which cannot overflow (w^2u^2/(1+w^2u^2) is
% written 1-1/q so it stays finite).

if nargin<7, idx=[]; end
if nargin<8 || isempty(nx), nx=0; end

wu=w*u;
q=1+wu.^2;
theta0=sum(k.*atan(wu)+l.*wu./q,1)/2;
log_rho=sum(k.*log(q),1)/4+sum(l.*(1-1./q),1)/2+(u*s).^2/8;
Z=exp(1i*theta0-log_rho);

if strcmpi(output,'cdf')
    H=-1i*Z./u;
elseif strcmpi(output,'pdf')
    H=Z;
else
    dx=(-(1i*u/2)).^nx;            % nx x-derivatives (R2)
    if strcmpi(output,'dens')
        H=dx.*Z;
        return
    elseif strcmpi(output,'k_deriv')
        ell=-0.5*log(1-1i*w(idx)*u);
    elseif strcmpi(output,'kk_deriv')
        ell=-0.5*log(1-1i*w(idx(1))*u).*(-0.5*log(1-1i*w(idx(2))*u));
    end
    H=1i*ell.*dx.*Z./u;
end
end
