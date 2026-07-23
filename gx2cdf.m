function [p,p_err,x_grid]=gx2cdf(x,w,k,lambda,s,m,varargin)

% GX2CDF Returns the cdf of a generalized chi-squared distribution.
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
% p=gx2cdf(x,w,k,lambda,s,m)
% p=gx2cdf(x,w,k,lambda,s,m,'upper')
% p=gx2cdf(x,w,k,lambda,s,m,'method','imhof','AbsTol',0,'RelTol',1e-7)
% [p,p_err]=gx2cdf(x,w,k,lambda,s,m,'method','ray','n_rays',1e7)
% [p,~,x_grid]=gx2cdf('full',w,k,lambda,s,m)
% etc.
%
% Example:
% p=gx2cdf(25,[1 -5 2],[1 2 3],[2 3 7],0,5)
%
% Required inputs:
% x         array of points at which to evaluate the cdf.
%           set x to 'full' to use IFFT method to return p over an array of
%           x that spans the distribution.
%           if method is 'ellipse' and 'x_scale' is 'log', these are log10
%           values of points measured from the finite tail, i.e.
%           log10(abs(x-m)), and will return log10 of p.
%
% w         row vector of weights of the non-central chi-squares
% k         row vector of degrees of freedom of the non-central chi-squares
% lambda    row vector of non-centrality paramaters (sum of squares of
%           means) of the non-central chi-squares
% s         scale of normal term
% m         offset
%
% Optional positional input:
% 'upper'   more accurate estimate of the complementary cdf when it's small.
%           if 'method' is 'tail', this option returns the upper infinite tail
%           approximation of the ccdf, otherwise the lower infinite tail approximation of the cdf.
%
% Optional name-value inputs:
% method    'auto' (default) tries to pick the best method for the parameters.
%           'imhof' for Imhof-Davies method.
%           'ray' for ray-trace method.
%           'ifft' for IFFT method.
%           'ruben' for Ruben's method. All w must be same sign and s=0.
%           'tail' for Das's infinite-tail approximation.
%           'pearson' for Imhof's extension to Pearson's 3-moment
%           approximation, extended again to include m and s.
%           'ellipse' for ellipse approximation. All w must be same sign and s=0.
% precision (for Imhof method) 'basic' (default) uses double precision, 'vpa' uses variable precision.
%           (for ray method) 'basic' uses double precision, 'log' (default)
%           returns the log10 of outputs smaller than realmin, 'vpa' uses variable precision.
% AbsTol    absolute error tolerance for the output. Default=1e-10.
% RelTol    relative error tolerance for the output. Default=1e-6.
%           AbsTol and RelTol are used only for Imhof's method, and ray method with grid integration.
%           The absolute OR the relative tolerance is satisfied.
%
% Options for ray method:
% force_mc  true to force Monte-Carlo instead of grid integration over
%           rays. Default=false.
% n_rays    no. of rays for Monte-Carlo integration. Larger=more accurate.
% gpu_batch no. of rays to compute on each GPU batch. 0 to use CPU.
%
% Options for Ruben method:
% n_ruben   no. of terms in Ruben's series. Larger=more accurate. Default=1000.
%
% Options for ifft method:
% span      span of x over which to compute the IFFT. Larger=more accurate.
% n_grid    number of grid points used for IFFT. Larger=more accurate. Default=1e6.
%
% Options for ellipse method:
% x_scale   'linear' (default). 'log' if input x is log10 values of x,
%           measured from the finite tail m of an elliptic distribution.
%           This will return log10 values of the probability.
%
% Outputs:
% p         computed cdf or ccdf.
%           (if ray method with 'precision' set to 'log') log10 of values < realmin
%           (if ray method with 'precision' set to 'vpa') symbolic output for values < realmin
%           (if tail method) cdf for lower tail, ccdf for upper tail. log10 of p for values < realmin.
%           
%           (if ellipse method with 'log' x_scale) log10 of p
% p_err     (if ray method with Monte-Carlo integration) standard error of the output p
%           (if Ruben's method) upper error bound of the output p
%           (if Imhof's method) logical array, true for outputs too close to 0 or 1 to compute exactly with
%           default settings. Try stricter tolerances.
%           (if ellipse method) relative error of p
%           (if ellipse method with 'log' x_scale) log10 of relative error of p
% x_grid    if input x is 'full', this returns the x points at which p is returned
%
% See also:
% <a href="matlab:open(strcat(fileparts(which('gx2cdf')),filesep,'doc',filesep,'GettingStarted.mlx'))">Getting Started guide</a>

parser = inputParser;
parser.KeepUnmatched = true;
addRequired(parser,'x',@(x) isreal(x) || strcmpi(x,'full'));
addRequired(parser,'w',@(x) isreal(x) && isrow(x));
addRequired(parser,'k',@(x) isreal(x) && isrow(x));
addRequired(parser,'lambda',@(x) isreal(x) && isrow(x));
addRequired(parser,'s',@(x) isreal(x) && isscalar(x));
addRequired(parser,'m',@(x) isreal(x) && isscalar(x));
addOptional(parser,'side','lower',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );
addParameter(parser,'method','auto');

parse(parser,x,w,k,lambda,s,m,varargin{:});
method=parser.Results.method;
side=parser.Results.side;

% initialize optional outputs so all dispatch paths assign them
p_err=[];
x_grid=[];

if strcmpi(x,'full')
    method='ifft';
end

if strcmpi(method,'auto')
    w_unique=unique(w);
    if ~s && isscalar(w_unique) % no s and only one unique weight
        % ncx2 fallback
        if (sign(w_unique)==1 && strcmpi(side,'lower')) || (sign(w_unique)==-1 && strcmpi(side,'upper'))
            p=ncx2cdf((x-m)/w_unique,sum(k),sum(lambda));
        else
            p=ncx2cdf((x-m)/w_unique,sum(k),sum(lambda),'upper');
        end
    elseif sum(abs(w))==0 && s % only normal term
        if strcmpi(side,'lower')
            p=normcdf(x,m,s);
        elseif strcmpi(side,'upper')
            p=normcdf(x,m,s,'upper');
        end
    elseif ~s
        if (all(w>0) && strcmpi(side,'lower'))||(all(w<0) && strcmpi(side,'upper'))  % no s and w same sign
            try
                [p,p_err]=gx2_ruben(x,w,k,lambda,m,varargin{:});
            catch
                [p,p_err]=gx2_imhof(x,w,k,lambda,0,m,varargin{:});
            end
        else
            [p,p_err]=gx2_imhof(x,w,k,lambda,s,m,varargin{:});
        end
    else
        [p,p_err]=gx2_imhof(x,w,k,lambda,s,m,varargin{:});
    end
elseif strcmpi(method,'ifft')
    [p,x_grid]=gx2_ifft(x,w,k,lambda,s,m,varargin{:},'output','cdf');
    p_err=[];
elseif strcmpi(method,'ray')
    [p,p_err]=gx2cdf_ray(x,w,k,lambda,s,m,varargin{:});
elseif strcmpi(method,'imhof')
    [p,p_err]=gx2_imhof(x,w,k,lambda,s,m,varargin{:});
elseif strcmpi(method,'ruben')
    gx2_check_method(method,w,s);
    [p,p_err]=gx2_ruben(x,w,k,lambda,m,varargin{:});
elseif strcmpi(method,'tail')
    p=gx2_tail(x,w,k,lambda,s,m,varargin{:});
elseif strcmpi(method,'pearson')
    p=gx2_pearson(x,w,k,lambda,s,m,varargin{:});
elseif strcmpi(method,'ellipse')
    gx2_check_method(method,w,s);
    [p,p_err]=gx2_ellipse(x,w,k,lambda,m,varargin{:});
end