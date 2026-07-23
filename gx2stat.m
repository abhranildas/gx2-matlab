function [mu,v,mode]=gx2stat(w,k,lambda,s,m)

    % GX2STAT Returns the mean, variance, and (optionally) mode of a generalized
    % chi-squared distribution.
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
    % [mu,v]=gx2stat(w,k,lambda,s,m)
    % [mu,v,mode]=gx2stat(w,k,lambda,s,m)
    %
    % Example:
    % [mu,v]=gx2stat([1 -5 2],[1 2 3],[2 3 7],0,5)
    % [mu,v,mode]=gx2stat([1 -5 2],[1 2 3],[2 3 7],0,5)
    %
    % Required inputs:
    % w         row vector of weights of the non-central chi-squares
    % k         row vector of degrees of freedom of the non-central chi-squares
    % lambda    row vector of non-centrality paramaters (sum of squares of
    %           means) of the non-central chi-squares
    % s         scale of normal term
    % m         offset
    %
    % Outputs:
    % mu        mean
    % v         variance
    % mode      mode. Unlike the mean and variance, it has no closed form, so
    %           it is located numerically as the root of the density derivative
    %           f'(x)=0 (using the analytic f' from gx2_dens_deriv), seeded at
    %           the mean. It is computed only when requested as the third
    %           output, so callers wanting just the mean/variance pay nothing.
    %           Assumes a single interior peak; for a density with no interior
    %           mode it falls back to maximizing the pdf on a wide interval.
    %
    % See also:
    % <a href="matlab:open(strcat(fileparts(which('gx2cdf')),filesep,'doc',filesep,'GettingStarted.mlx'))">Getting Started guide</a>

    parser = inputParser;
    addRequired(parser,'w',@(x) isreal(x) && isrow(x));
    addRequired(parser,'k',@(x) isreal(x) && isrow(x));
    addRequired(parser,'lambda',@(x) isreal(x) && isrow(x));
    addRequired(parser,'s',@(x) isreal(x) && isscalar(x));
    addRequired(parser,'m',@(x) isreal(x) && isscalar(x));
    parse(parser,w,k,lambda,s,m);

    mu=dot(w,k+lambda)+m;
    v=2*dot(w.^2,k+2*lambda)+s^2;

    if nargout>2
        % Mode: root of f'(x)=0, seeded at the mean. fzero brackets a sign
        % change of f' on its own; the density derivative comes from
        % gx2_dens_deriv. Tiny tail-density evaluations can clip (a benign
        % gx2_imhof warning) -- silence just that id during the search. If no
        % interior stationary point is found (a monotone density, mode at the
        % support edge), fall back to maximizing the pdf on a wide interval.
        wstate=warning('off','gx2:imhofClip'); cln=onCleanup(@() warning(wstate)); %#ok<NASGU>
        fp=@(xx) gx2_dens_deriv(xx,w,k,lambda,s,m,1);
        try mode=fzero(fp,mu); catch, mode=NaN; end
        if ~isfinite(mode)
            sd=sqrt(v);
            mode=fminbnd(@(xx) -gx2pdf(xx,w,k,lambda,s,m),mu-8*sd,mu+8*sd);
        end
    end