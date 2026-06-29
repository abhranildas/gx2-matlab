function [mu,v]=gx2stat(w,k,lambda,s,m)

    % GX2STAT Returns the mean and variance of a generalized chi-squared distribution.
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
    % [mu,v]=gx2stat(w,k,lambda,s,m)
    %
    % Example:
    % [mu,v]=gx2stat([1 -5 2],[1 2 3],[2 3 7],0,5)
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