function [grad,hess]=norm_err_grad_bd(mu0,v0,mu1,v1,quad,varargin)

% NORM_ERR_GRAD_BD Gradient (and Hessian) of the total two-class
% classification error with respect to the boundary coefficients, for a
% quadratic boundary shared between two normal classes.
%
% The error is E = p0*P(q(x)<=0 | x~N(mu0,v0)) + p1*P(q(x)>0 | x~N(mu1,v1)),
% i.e. class 0 mass on the class-1 side of the boundary plus class 1 mass on
% the class-0 side, with the convention that class 0 is favored where
% q(x)>0 (matching IntClassNorm's norm_class_opt_bd). Since
% P(q(x)>0)=1-F(0), this is a signed combination of the two classes'
% CDF_GRAD_BD outputs at x0=0: dE = p0*dF0 - p1*dF1 for the gradient, and
% likewise block by block for the Hessian.
%
% If the shared boundary lands exactly on one class's density cusp (see
% CDF_GRAD_BD), it generically lands on both classes' cusps at once
% (the two F's share the same boundary coefficients), and resolving each
% class's divergence independently before subtracting can turn a genuine
% inf-inf cancellation into nan -- the two classes' divergences may be
% mirror images that cancel in the combined error even though neither
% class's own Hessian is finite there. So this cusp is resolved on the
% *combined* quantity directly: the same shrinking-probe test as
% CDF_GRAD_BD, applied to p0*F0-p1*F1 as a whole rather than to each
% class alone.
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
% grad=norm_err_grad_bd(mu0,v0,mu1,v1,quad)
% [grad,hess]=norm_err_grad_bd(mu0,v0,mu1,v1,quad,'p0',0.3,'p1',0.7)
%
% Required inputs:
% mu0,v0    mean (column vector) and covariance of class 0
% mu1,v1    mean (column vector) and covariance of class 1
% quad      struct with the boundary's quadratic form coefficients, shared
%           between the two classes (e.g. the optimal quadratic
%           discriminant between them) -- see CDF_GRAD_BD
%
% Optional name-value inputs:
% p0,p1     class priors (weights on each class's error term). Default 0.5 each.
% wrt, AbsTol, RelTol, precision, n_ruben
%           forwarded to CDF_GRAD_BD's underlying computation for
%           each class -- see there for defaults and meaning.
%
% Outputs:
% grad      gradient of E, same fields as CDF_GRAD_BD's gradient
% hess      Hessian of E, only if a second output is requested -- same
%           fields as CDF_GRAD_BD's Hessian
%
% See also:
% cdf_grad_bd

parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'mu0',@isnumeric);
addRequired(parser,'v0',@isnumeric);
addRequired(parser,'mu1',@isnumeric);
addRequired(parser,'v1',@isnumeric);
addRequired(parser,'quad',@isstruct);
groups={'q2','q1','q0'};
addParameter(parser,'p0',0.5,@(x) isreal(x) && isscalar(x));
addParameter(parser,'p1',0.5,@(x) isreal(x) && isscalar(x));
addParameter(parser,'wrt',groups,@(c) iscell(c) && all(ismember(lower(c),groups)));
addParameter(parser,'AbsTol',1e-10,@(x) isreal(x) && isscalar(x) && (x>=0));
addParameter(parser,'RelTol',1e-6,@(x) isreal(x) && isscalar(x) && (x>=0));
addParameter(parser,'precision','basic',@(x) strcmpi(x,'basic')||strcmpi(x,'vpa'));
addParameter(parser,'n_ruben',1e3,@(x) isscalar(x) && (x>0) && (x==round(x)));
parse(parser,mu0,v0,mu1,v1,quad,varargin{:});

p0=parser.Results.p0;
p1=parser.Results.p1;
wrt=parser.Results.wrt;
AbsTol=parser.Results.AbsTol;
RelTol=parser.Results.RelTol;
precision=parser.Results.precision;
n_ruben=parser.Results.n_ruben;

want_hess=nargout>=2;
[g0,h0,at_cusp0,scale0,inv0,nx]=boundary_raw(0,mu0,v0,quad,wrt,want_hess,AbsTol,RelTol,precision,n_ruben);
[g1,h1,at_cusp1,scale1,inv1,~ ]=boundary_raw(0,mu1,v1,quad,wrt,want_hess,AbsTol,RelTol,precision,n_ruben);

raw_grad=struct();
gfields=fieldnames(g1);
for fi=1:numel(gfields)
    ky=gfields{fi};
    raw_grad.(ky)=p0*g0.(ky)-p1*g1.(ky);
end

raw_hess=[];
if want_hess
    raw_hess=struct();
    hfields=fieldnames(h1);
    for fi=1:numel(hfields)
        ky=hfields{fi};
        raw_hess.(ky)=p0*h0.(ky)-p1*h1.(ky);
    end
end

at_cusp=at_cusp0|at_cusp1;
if any(at_cusp)
    cusp_scale=max(scale0,scale1);
    inversion_grad_hess=@(s_eff,want_hess_probe) combined_probe(s_eff,want_hess_probe,inv0,inv1,p0,p1);
    [raw_grad,raw_hess]=resolve_cusp_terms(raw_grad,raw_hess,at_cusp,inversion_grad_hess,cusp_scale,nx);
end

[grad,hess]=finalize_grad_hess(raw_grad,raw_hess,want_hess,nx);

end

function [cg,ch]=combined_probe(s_eff,want_hess_probe,inv0,inv1,p0,p1)
% p0*F0-p1*F1 at an artificial normal term s_eff, applied to both classes'
% own inversion_grad_hess probes -- see NORM_ERR_GRAD_BD's header for why
% this is combined *before* resolving rather than after.
[cg0,ch0]=inv0(s_eff,want_hess_probe);
[cg1,ch1]=inv1(s_eff,want_hess_probe);
cg=struct();
gfields=fieldnames(cg1);
for fi=1:numel(gfields)
    ky=gfields{fi};
    cg.(ky)=p0*cg0.(ky)-p1*cg1.(ky);
end
ch=[];
if want_hess_probe
    ch=struct();
    hfields=fieldnames(ch1);
    for fi=1:numel(hfields)
        ky=hfields{fi};
        ch.(ky)=p0*ch0.(ky)-p1*ch1.(ky);
    end
end
end
