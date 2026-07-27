function [p,p_err]=gx2_ruben(x,w,k,l,m,varargin)

% GX2_RUBEN Ruben's series method for the generalized chi-square cdf/pdf.
% Requires all w the same sign.
%
% Abhranil Das
% Center for Perceptual Systems, University of Texas at Austin
% Comments, questions, bugs to abhranil.das@utexas.edu
%
% A thin wrapper over gx2_ruben_coeffs (the x-independent series setup) and
% gx2_ruben_eval (evaluation at given points); this function's public
% behavior is unchanged from before the split. Callers who evaluate many
% points against the same (w,k,l) -- e.g. gx2_dens_deriv's mixed-sign
% convolution, sampling one scalar point per quadrature node -- should call
% gx2_ruben_coeffs once and reuse it via gx2_ruben_eval directly, rather than
% going through this wrapper (which rebuilds the coefficients every call).
%
% Usage:
% p=gx2_ruben(x,w,k,l,m)
% p=gx2_ruben(x,w,k,l,m,'upper')
% p=gx2_ruben(x,w,k,l,m,'output','pdf','nx',1)
%
% Required inputs:
% x         array of points to evaluate at
% w         row vector of weights, all the same sign
% k         row vector of degrees of freedom
% l         row vector of non-centralities
% m         scalar offset
%
% Optional positional/name-value inputs:
% side      'lower' (default) or 'upper'
% output    'cdf' (default) or 'pdf'
% n_ruben   cap on the number of series terms. Default=1e3.
% nx        x-derivative order (pdf output only). Default=0.
%
% Outputs:
% p         cdf or pdf values, same size as x
% p_err     truncation-error indicator
%
% See also:
% gx2_ruben_coeffs, gx2_ruben_eval

parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'x',@(x) isreal(x));
addRequired(parser,'w',@(x) isreal(x) && isrow(x)  && (all(x>0)||all(x<0)) );
addRequired(parser,'k',@(x) isreal(x) && isrow(x));
addRequired(parser,'l',@(x) isreal(x) && isrow(x));
addRequired(parser,'m',@(x) isreal(x) && isscalar(x));
addOptional(parser,'side','lower',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );
addParameter(parser,'output','cdf',@(x) strcmpi(x,'cdf') || strcmpi(x,'pdf') );
addParameter(parser,'n_ruben',1e3,@(x) isscalar(x) && (x>0) && (x==round(x)));
addParameter(parser,'nx',0,@(x) isscalar(x) && (x>=0) && (x==round(x)));

parse(parser,x,w,k,l,m,varargin{:});
side=parser.Results.side;
output=parser.Results.output;
n_ruben=parser.Results.n_ruben;
nx=parser.Results.nx;

coeffs=gx2_ruben_coeffs(w,k,l,'n_ruben',n_ruben);
[p,p_err]=gx2_ruben_eval(coeffs,x,m,side,'output',output,'nx',nx);

end
