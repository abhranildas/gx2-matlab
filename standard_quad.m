function quad_s=standard_quad(quad,mu,v)
% standardize quadratic coefficients
sv=sqrtm(v); % symmetric square root, computed once
quad_s.q2=sv*quad.q2*sv;
quad_s.q1=sv*(2*quad.q2*mu+quad.q1);
quad_s.q0=mu'*quad.q2*mu+quad.q1'*mu+quad.q0;