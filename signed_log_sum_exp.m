function s = signed_log_sum_exp(logs, dim)
% SIGNED_LOG_SUM_EXP   Log-space summation with signed inputs.
%
%   s = SIGNED_LOG_SUM_EXP(logs,dim) returns the signed log–sum–exp of the numbers
%   whose values are represented by 'logs' using the following convention:
%       For a positive number x, represent it as  -log10(x) (which is negative).
%       For a negative number x, represent it as  log10(|x|) (which is positive).
%   The output s uses the same convention:
%       s < 0  <==>  S = 10^(-|s|) > 0,
%       s > 0  <==>  S = -10^(-|s|) < 0.
%
%   (If r==0 (exact cancellation) you may wish to handle it specially.)
%
%   Example:
%       % To compute  10^-3 - 10^-4, represent 10^-3 (positive) as -3 
%       % and 10^-4 (to be subtracted) as +4.
%       logs = [-3, 4];
%       s = signed_log_sum_exp(logs)
%       % Here, the sum is 10^-3 - 10^-4 = 0.9e-3, whose representation is:
%       %   -log10(0.9e-3) = -( -3 + log10(0.9) ) ≈ -3.0458.


if nargin < 2
    dim = 'all';
end

if isempty(logs)
    s = -inf;
    return;
end

logs(isnan(logs))=-inf;

% Separate magnitude and sign.
mag = abs(logs);
term_sign = -sign(logs);   % positive for logs<0, negative for logs>0

% Factor out the dominant (largest) term.
m = min(mag, [], dim);

% Compute the rescaled sum:
r = sum(term_sign .* 10.^ ( -(mag - m) ), dim);

% (Optional) Handle complete cancellation:
if any(r == 0)
    % You might choose to return a special value (for example, NaN)
    r(r==0) = eps; % or, alternatively, s(r==0)=NaN and exit.
end

% Reconstruct the signed log of the sum.
s = zeros(size(r));
% If r>0, then S>0 and the signed log must be negative:
pos = r > 0;
s(pos) = -( m(pos) - log10( r(pos) ) );
% If r<0, then S<0 and the signed log must be positive:
neg = r < 0;
s(neg) = m(neg) - log10( abs(r(neg)) );

s(m==inf)=-inf;
end
