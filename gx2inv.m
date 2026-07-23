function x=gx2inv(p,w,k,lambda,s,m,varargin)

    % GX2INV Returns the inverse cdf of a generalized chi-squared distribution.
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
    addRequired(parser,'w',@(x) isreal(x) && isrow(x));
    addRequired(parser,'k',@(x) isreal(x) && isrow(x));
    addRequired(parser,'lambda',@(x) isreal(x) && isrow(x));
    addRequired(parser,'s',@(x) isreal(x) && isscalar(x));
    addRequired(parser,'m',@(x) isreal(x) && isscalar(x));
    addOptional(parser,'side','lower',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );

    parse(parser,p,w,k,lambda,s,m,varargin{:});

    side=parser.Results.side;

    w_unique=unique(w);
    if ~s && isscalar(w_unique) && all(p>0)
        % native ncx2 fallback
        if strcmpi(side,'upper')
            p=1-p;
        end
        if sign(w_unique)==1
            x=ncx2inv(p,sum(k),sum(lambda))*w_unique+m;
        elseif sign(w_unique)==-1
            x=ncx2inv(1-p,sum(k),sum(lambda))*w_unique+m;
        elseif w_unique==0
            x=0;
        end
    else
        [mu,v]=gx2stat(w,k,lambda,s,m);
        sd=sqrt(v);
        if all(p>0)
            % Safeguarded Newton on gx2cdf, using the analytic pdf as the
            % derivative (dG/dx = +-f). Bracketed with a bisection fallback, so
            % it is as robust as fzero but needs fewer gx2cdf evaluations. All
            % requested probabilities are solved simultaneously: each Newton
            % step evaluates the cdf and pdf at the whole vector of iterates in
            % one shared Imhof quadrature sweep, so the integral cost scales
            % with the number of iterations, not with numel(p). The pdf takes
            % the same passthrough options as the cdf, minus the 'upper'/'lower'
            % side (captured separately); the derivative sign is -f for the
            % upper tail. The log-probability branch stays on fzero.
            pdfargs=namedargs2cell(parser.Unmatched);
            dsgn=1; if strcmpi(side,'upper'), dsgn=-1; end
            G =@(xx) gx2cdf(xx,w,k,lambda,s,m,varargin{:});
            dG=@(xx) dsgn*gx2pdf(xx,w,k,lambda,s,m,pdfargs{:});
            x=gx2inv_newton(G,dG,p,mu,sd);
        else % log probability, inverted using sym
            x=arrayfun(@(p) fzero(@(x) log_gx2cdf(x,w,k,lambda,s,m,varargin{:})-p,mu),p);
        end
    end

function x=gx2inv_newton(G,dG,p,x0,sd)
    % Roots of the monotone residual r(x)=G(x)-p by Newton with bisection
    % fallback (the classic rtsafe scheme), solving every element of p at once.
    % G and dG accept a vector of iterates and return the cdf/pdf at all of them
    % in one Imhof sweep, so each iteration costs two integrals for the whole
    % batch rather than two per probability. xl/xh label the bracket ends by the
    % sign of r, so it works whether G increases (lower tail) or decreases
    % (upper tail). Any element that cannot be bracketed falls back to fzero.
    % Intermediate cdf/pdf evaluations can clip tiny tail values (a benign
    % gx2_imhof warning); silence just that id during the solve -- the returned
    % roots are still governed by convergence, and extreme tails use the fzero
    % log-probability branch instead. Auto-restored on return.
    wstate=warning('off','gx2:imhofClip'); cln=onCleanup(@() warning(wstate));
    tol=4*eps; maxit=100;
    if ~isfinite(sd) || sd<=0, sd=1; end

    % Bracket every probability with one symmetric window about x0, widened
    % until a sign change straddles the residual for all of them. The ends a,b
    % stay scalar during expansion (fa,fb are shaped like p); monotonicity of G
    % keeps a bracket valid once found, so widening never breaks earlier ones.
    step=max(sd,1e-3); a=x0-step; b=x0+step; fa=G(a)-p; fb=G(b)-p; it=0;
    while any(fa(:).*fb(:)>0) && it<80
        step=1.6*step; a=a-step; b=b+step; fa=G(a)-p; fb=G(b)-p; it=it+1;
    end
    unbr=fa.*fb>0;                               % elements never bracketed
    lowa=fa<0;                                   % where a is the low end of r
    xl=a*lowa+b*(~lowa); xh=b*lowa+a*(~lowa);

    x=0.5*(a+b)*ones(size(p)); dxold=abs(b-a)*ones(size(p)); dx=dxold;
    fx=G(x)-p; dfx=dG(x); act=~unbr;             % iterate the bracketed ones
    for it=1:maxit
        % per-element choice: bisect when Newton leaves the bracket, stalls, or
        % has a bad derivative; otherwise take the Newton step
        bis=((x-xh).*dfx-fx).*((x-xl).*dfx-fx)>0 | abs(2*fx)>abs(dxold.*dfx) ...
            | dfx==0 | ~isfinite(dfx);
        dxold=dx;
        dxb=0.5*(xh-xl); xb=xl+dxb;              % bisection candidate
        dxn=fx./dfx;     xn=x-dxn;               % Newton candidate
        dx=dxn; dx(bis)=dxb(bis);                % select (indexing avoids 0*NaN)
        x =xn;  x(bis) =xb(bis);
        if all(abs(dx(act))<=tol*(1+abs(x(act)))), break; end
        fx=G(x)-p; dfx=dG(x);
        lo=fx<0; xl(lo)=x(lo); xh(~lo)=x(~lo);   % tighten the bracket
    end

    bad=unbr | ~isfinite(x);                     % scalar fzero for stragglers
    for j=find(bad(:))'
        x(j)=fzero(@(xx) G(xx)-p(j),x0);
    end