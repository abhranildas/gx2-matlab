function phi=gx2char(t,w,k,lambda,s,m)

    % GX2CHAR Returns the characteristic function of a generalized chi-squared distribution.
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
    % phi=gx2char(t,w,k,lambda,s,m)
    %
    % Example:
    % phi=gx2char(linspace(-10,10,100),[1 -5 2],[1 2 3],[2 3 7],0,5)
    %
    % Required inputs:
    % t         array of points at which to compute the characteristic function
    % w         row vector of weights of the non-central chi-squares
    % k         row vector of degrees of freedom of the non-central chi-squares
    % lambda    row vector of non-centrality paramaters (sum of squares of
    %           means) of the non-central chi-squares
    % s         scale of normal term
    % m         offset
    %
    % Outputs:
    % phi       characteristic function
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

    t_flat=t(:); % flatten input array

    phi=exp(1i*m*t_flat+1i*t_flat.*sum((w.*lambda)./(1-2i*t_flat.*w),2)-s^2*t_flat.^2/2)./...
        prod((1-2i*w.*t_flat).^(k/2),2);

    % reshape output to shape of input t
    phi=reshape(phi,size(t));

