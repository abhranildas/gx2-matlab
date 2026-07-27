function coeffs=gx2_ruben_coeffs(w,k,l,varargin)

% GX2_RUBEN_COEFFS Ruben's series expansion coefficients for the generalized
% chi-square with weights w, degrees of freedom k and non-centralities l.
%
% Abhranil Das
% Center for Perceptual Systems, University of Texas at Austin
% Comments, questions, bugs to abhranil.das@utexas.edu
%
% Depends only on (w,k,l) and the term-count cap n_ruben -- not on the
% evaluation point x or offset m -- so callers evaluating many points against
% the same (w,k,l) (e.g. gx2_dens_deriv's mixed-sign convolution, which
% samples one scalar point per quadrature callback) should compute this once
% and reuse it via gx2_ruben_eval, rather than recomputing the series from
% scratch at every point.
%
% Usage:
% coeffs=gx2_ruben_coeffs(w,k,l)
% coeffs=gx2_ruben_coeffs(w,k,l,'n_ruben',500)
%
% Required inputs:
% w         row vector of weights, all the same sign
% k         row vector of degrees of freedom
% l         row vector of non-centralities
%
% Optional name-value inputs:
% n_ruben   cap on the number of series terms. Default=1e3.
%
% Output:
% coeffs    struct (fields a, N, beta, M, w_pos) for use with gx2_ruben_eval
%
% See also:
% gx2_ruben_eval, gx2_ruben

parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'w',@(x) isreal(x) && isrow(x)  && (all(x>0)||all(x<0)) );
addRequired(parser,'k',@(x) isreal(x) && isrow(x));
addRequired(parser,'l',@(x) isreal(x) && isrow(x));
addParameter(parser,'n_ruben',1e3,@(x) isscalar(x) && (x>0) && (x==round(x)));

parse(parser,w,k,l,varargin{:});
n_ruben=parser.Results.n_ruben;

w_pos=true;
if all(w<0)
    w=-w; w_pos=false;
end
beta=0.90625*min(w);
M=sum(k);
n=(1:n_ruben-1)';

% compute the g's
g=sum(k.*(1-beta./w).^n,2)+ beta*n.*((1-beta./w).^(n-1))*(l./w)';

% compute the expansion coefficients, stopping once the leftover series mass
% is negligible. The a_j are nonnegative and sum to 1, so the tail mass
% 1-sum(a_{1:N}) both bounds the truncation error and decreases monotonically.
% The stop uses only this cheap coefficient recurrence -- not the chi-square
% grid evaluated by gx2_ruben_eval -- so the term count N is fixed in a single
% pass, and the expensive evaluation is then done once at that reduced size.
% n_ruben is the safety cap; most cases converge in ~10^2 terms well under it.
masstol=1e-14;
a=nan(n_ruben,1);
a(1)=sqrt(exp(-sum(l))*beta^M*prod(w.^(-k)));
if a(1)<realmin
    error('Underflow error: some series coefficients are smaller than machine precision.')
end
cum=a(1); N=n_ruben;
for j=1:n_ruben-1
    a(j+1)=dot(flip(g(1:j)),a(1:j))/(2*j);
    cum=cum+a(j+1);
    if 1-cum<masstol, N=j+1; break; end
end
a=a(1:N);

coeffs=struct('a',a,'N',N,'beta',beta,'M',M,'w_pos',w_pos);

end
