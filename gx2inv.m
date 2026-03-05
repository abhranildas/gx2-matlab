function x=gx2inv(p,w,k,lambda,s,m,varargin)

    % GX2INV Returns the inverse cdf of a generalized chi-squared distribution.
    %
    % Abhranil Das
    % Center for Perceptual Systems, University of Texas at Austin
    % Comments, questions, bugs to abhranil.das@utexas.edu
    % If you use this code, please cite:
    % 1. <a href="matlab:web('https://arxiv.org/abs/2012.14331')"
    % >A method to integrate and classify normal distributions</a>
    % 2. <a href="matlab:web('https://arxiv.org/abs/2404.05062')"
    % >New methods for computing the generalized chi-square distribution</a>
    %
    % Usage:
    % x=gx2inv(p,w,k,lambda,s,m)
    % x=gx2inv(p,w,k,lambda,s,m,'upper','method','imhof','AbsTol',0,'RelTol',1e-7)
    % etc.
    %
    % Example:
    % x=gx2inv(0.9,[1 -5 2],[1 2 3],[2 3 7],0,5)
    % x=gx2inv(-100,[1 -5 2],[1 2 3],[2 3 7],0,5,'upper','method','ray','n_rays',1e4)
    %
    % Required inputs:
    % p         probabilities at which to evaluate the inverse cdf.
    %           Negative values indicate log probability, that can be used
    %           to invert probabilities < realmin, using ray, ellipse, or tail cdf methods.
    % w         row vector of weights of the non-central chi-squares
    % k         row vector of degrees of freedom of the non-central chi-squares
    % lambda    row vector of non-centrality paramaters (sum of squares of
    %           means) of the non-central chi-squares
    % s         scale of normal term
    % m         offset
    %
    % Optional positional input:
    % 'upper'   for more accurate quantiles when entering an upper tail
    %           probability (complementary cdf)
    %
    % Optional name-value inputs:
    % This function numerically finds roots of the gx2cdf function, so most
    % options for the gx2cdf function can be used here, eg 'method' and
    % 'x_scale', which will be passed on to gx2cdf
    %
    % Output:
    % x         computed quantile
    %
    % See also:
    % <a href="matlab:open(strcat(fileparts(which('gx2cdf')),filesep,'doc',filesep,'GettingStarted.mlx'))">Getting Started guide</a>

    parser=inputParser;
    parser.KeepUnmatched=true;
    addRequired(parser,'p',@(x) isreal(x) && all(x<=1));
    addRequired(parser,'w',@(x) isreal(x));
    addRequired(parser,'k',@(x) isreal(x));
    addRequired(parser,'lambda',@(x) isreal(x));
    addRequired(parser,'s',@(x) isreal(x) && isscalar(x));
    addRequired(parser,'m',@(x) isreal(x) && isscalar(x));
    addOptional(parser,'side','lower',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );

    parse(parser,p,w,k,lambda,s,m,varargin{:});

    side=parser.Results.side;

    if ~s && isscalar(unique(w)) && all(p>0)
        % native ncx2 fallback
        if strcmpi(side,'upper')
            p=1-p;
        end
        if sign(unique(w))==1
            x=ncx2inv(p,sum(k),sum(lambda))*unique(w)+m;
        elseif sign(unique(w))==-1
            x=ncx2inv(1-p,sum(k),sum(lambda))*unique(w)+m;
        elseif unique(w)==0
            x=0;
        end
    else
        mu=gx2stat(w,k,lambda,s,m);
        if p>0
            x=arrayfun(@(p) fzero(@(x) gx2cdf(x,w,k,lambda,s,m,varargin{:})-p,mu),p);
        else % log probability, inverted using sym
            x=arrayfun(@(p) fzero(@(x) log_gx2cdf(x,w,k,lambda,s,m,varargin{:})-p,mu),p);
        end
    end