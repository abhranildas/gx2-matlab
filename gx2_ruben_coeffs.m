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
% n_ruben   cap on the number of series terms. Default=1e3. If the series
%           hasn't converged within this many terms (checked cheaply from
%           the coefficients alone, before the expensive per-x evaluation
%           done by gx2_ruben_eval) -- which happens when the smallest |w|
%           is small enough to need far more terms than are efficient to
%           compute here -- this returns coeffs.a=NaN rather than a value
%           truncated at an arbitrary, unverified point. Callers that
%           dispatch across methods (gx2cdf's method='auto', and
%           gx2_dens_deriv) already treat a non-finite Ruben output as "fall
%           back to Imhof".
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
%
% Either branch below can instead exhaust n_ruben without its stopping
% criterion ever firing -- a real outcome, not a formal edge case: it happens
% whenever the smallest |w_j| is small enough that the implied scale
% beta=0.90625*min|w_j| makes (x-m)/beta huge, which is exactly what a
% near-rank-deficient quadratic boundary produces just off its exact
% zero-eigenvalue point. There the series provably still converges, but needs
% many more terms than are safe to compute here (a single extra order of
% magnitude in N costs two more in runtime, since the recursion is O(N^2));
% confirmed directly on such a case: it doesn't converge or match Imhof until
% N~2e4, taking >0.5s, versus Imhof's ~0.01s for the same point at any N. So
% this cap is a genuine efficiency boundary, not just a safety net, and a
% fixed larger cap would only move the failure to a slightly more extreme
% case, not remove it. We signal this with NaN (checked via the `converged`
% flag, since MATLAB has no for/else) rather than returning whatever
% badly-truncated partial sum the loop stopped at -- both of this function's
% callers already treat a non-finite Ruben output as "fall back to Imhof"
% (gx2cdf's auto method, and gx2_dens_deriv's same-sign and mixed-sign
% routes), so this reuses that existing, already-tested path instead of
% adding a new failure signal each caller would need to learn about
% separately.
masstol=1e-14;
a=nan(n_ruben,1);
a(1)=sqrt(exp(-sum(l))*beta^M*prod(w.^(-k)));
if a(1)<realmin
    % The true leading coefficient underflows (e.g. when some non-centrality
    % in l is large, as happens for a quadratic form whose curvature is small
    % relative to its linear part). The a_j are nonnegative and sum to 1, but
    % that overall scale is lost here -- only their relative sizes survive,
    % since the recursion below is linear and homogeneous in a(1:j). Recover
    % the coefficients up to that lost scale (starting from b(1)=1 instead of
    % the unrepresentable true a(1)), then renormalize at the end so they sum
    % to 1, exactly as they must.
    b=nan(n_ruben,1);
    b(1)=1;
    cum=b(1); N=n_ruben; converged=false;
    for j=1:n_ruben-1
        b(j+1)=dot(flip(g(1:j)),b(1:j))/(2*j);
        cum=cum+b(j+1);
        if b(j+1)<masstol*cum, N=j+1; converged=true; break; end
    end
    if ~converged
        coeffs=struct('a',nan,'N',1,'beta',beta,'M',M,'w_pos',w_pos);
        return
    end
    a=b(1:N)/cum;
else
    cum=a(1); N=n_ruben; converged=false;
    for j=1:n_ruben-1
        a(j+1)=dot(flip(g(1:j)),a(1:j))/(2*j);
        cum=cum+a(j+1);
        if 1-cum<masstol, N=j+1; converged=true; break; end
    end
    if ~converged
        coeffs=struct('a',nan,'N',1,'beta',beta,'M',M,'w_pos',w_pos);
        return
    end
    a=a(1:N);
end

coeffs=struct('a',a,'N',N,'beta',beta,'M',M,'w_pos',w_pos);

end
