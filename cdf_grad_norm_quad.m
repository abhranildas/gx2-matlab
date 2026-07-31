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
% n_ruben   term-count cap passed through to the Ruben-series density
%           derivatives used on the s==0 route (see gx2_dens_deriv, gx2_ruben).
%           Only affects that route; ignored when s~=0. Default=1e3.
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
% hess      struct, only if a second output is requested: the six blocks
%           q0q0, q0q1, q0q2, q1q1, q1q2 (a 3-tensor) and q2q2 (a 4-tensor),
%           each carrying a trailing dimension per threshold. If the shared
%           boundary lands exactly on this class's density cusp (an
%           antipodal-mean, swapped-covariance-type coincidence -- see
%           gx2_derivatives.md open item 3.3), the affected entries are
%           resolved to a correctly-signed inf rather than an arbitrary
%           finite number or nan.
%
% See also:
% cdf_grad_gx2, norm_err_grad_bd, norm_quad_to_gx2_params, gx2cdf, gx2pdf, gx2char

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
addParameter(parser,'n_ruben',1e3,@(x) isscalar(x) && (x>0) && (x==round(x)));
parse(parser,x,mu,v,quad,varargin{:});

wrt=parser.Results.wrt;
AbsTol=parser.Results.AbsTol;
RelTol=parser.Results.RelTol;
precision=parser.Results.precision;
n_ruben=parser.Results.n_ruben;

want_hess=nargout>=2;
[raw_grad,raw_hess,at_cusp,cusp_scale,inversion_grad_hess,nx]= ...
    boundary_raw(x,mu,v,quad,wrt,want_hess,AbsTol,RelTol,precision,n_ruben);

% Override only the exact-cusp points, via an independent always-valid
% method (see resolve_cusp_terms) -- every other point, and every other
% problem that never hits this measure-adjacent exact coincidence, is
% completely untouched.
if any(at_cusp)
    [raw_grad,raw_hess]=resolve_cusp_terms(raw_grad,raw_hess,at_cusp,inversion_grad_hess,cusp_scale,nx);
end

[grad,hess]=finalize_grad_hess(raw_grad,raw_hess,want_hess,nx);
