function [w,k,lambda,s,m,aux]=norm_quad_to_gx2_params(mu,v,quad,varargin)

% NORM_QUAD_TO_GX2_PARAMS A quadratic form of a normal vector is distributed
% as a generalized chi-squared. This function takes the multinormal parameters
% and the quadratic coefficients and returns the parameters of the generalized
% chi-squared.
%
% Abhranil Das
% Center for Perceptual Systems, University of Texas at Austin
% Comments, questions, bugs to abhranil.das@utexas.edu
% Thanks to bugfix by Florian Roth, florian.roth2@tu-dresden.de.
% If you use this code, please cite:
% 1. <a href="matlab:web('https://doi.org/10.1167/jov.21.10.1','-browser')"
% >A method to integrate and classify normal distributions</a>
% 2. <a href="matlab:web('https://www.tandfonline.com/doi/abs/10.1080/00949655.2025.2501401','-browser')"
% >New methods to compute the generalized chi-square distribution</a>
%
% Example:
% mu=[1;2]; % mean
% v=[2 1; 1 3]; % covariance matrix
% % Say q(x)=(x1+x2)^2-x1-1 = [x1;x2]'*[1 1; 1 1]*[x1;x2] + [-1;0]'*[x1;x2] - 1:
% quad.q2=[1 1; 1 1];
% quad.q1=[-1;0];
% quad.q0=-1;
%
% [w,k,lambda,s,m]=norm_quad_to_gx2_params(mu,v,quad)
%
% Required inputs:
% mu        column vector of normal mean
% v         normal covariance matrix
% quad      struct with quadratic form coefficients:
%               q2      matrix of quadratic coefficients
%               q1      column vector of linear coefficients
%               q0      scalar constant
%
% Optional name-value inputs:
% merge     true by default. Merges the non-central chi-square components
%           with close-enough weights (using uniquetol) into single
%           components. Set false to return all raw exact components (more
%           precise).
% Outputs:
% w         row vector of weights of the non-central chi-squares
% k         row vector of degrees of freedom of the non-central chi-squares
% lambda    row vector of non-centrality paramaters (sum of squares of
%           means) of the non-central chi-squares
% s         scale of normal term
% m         offset
% aux       (optional 6th output) struct with the eigen-structure of the
%           standardized quadratic S*Q2*S, reused by cdf_grad_norm_quad:
%           S (=Sigma^{1/2}), V (eigenvectors), d (eigenvalues, row), and
%           b (=V'*S*(2*Q2*mu+q1)). These are the full, pre-merge quantities.
%
% See also:
% <a href="matlab:open(strcat(fileparts(which('gx2cdf')),filesep,'doc',filesep,'GettingStarted.mlx'))">Getting Started guide</a>

% standardize the space:

parser = inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'mu',@isnumeric);
addRequired(parser,'v',@isnumeric);
addRequired(parser,'quad');
addParameter(parser,'merge',true,@islogical);

parse(parser,mu,v,quad,varargin{:});
merge=parser.Results.merge;

q2_sym=0.5*(quad.q2+quad.q2'); % symmetrize q2

% compute sqrtm of v while avoiding small negative eigenvalues
[R, D] = eig(v);
d = diag(D);
d(d < 0) = 0;              % Threshold any small negatives to zero
sqrt_d = sqrt(d);          % Compute the square roots
sqrt_v = R * diag(sqrt_d) * R';  % Reassemble the square root matrix

q2=sqrt_v*q2_sym*sqrt_v;
q2=(q2+q2')/2; % symmetrize q2 again
q1=sqrt_v*(2*q2_sym*mu+quad.q1);
q0=mu'*q2_sym*mu+quad.q1'*mu+quad.q0;

[R,D]=eig(q2);
d=diag(D)';
b=(R'*q1)';

% Split the eigenvalues into nonzero (chi-square directions) and effectively
% zero (purely linear directions, which contribute to the normal term s). This
% uses a relative rank tolerance rather than an exact d==0 test: a numerically
% tiny eigenvalue from a rank-deficient q2 would otherwise be kept as a
% chi-square with a near-zero weight and a blown-up non-centrality
% b^2/(4w^2), which overflows the Ruben series downstream.
dtol=max(abs(d))*numel(d)*eps;
nz=abs(d)>dtol;                       % nonzero (chi-square) eigenvalues

if merge
    [w,~,ic]=uniquetol(d(nz)); % unique non-zero eigenvalues
    k=accumarray(ic,1)'; % total dof of each eigenvalue

    % lambda=arrayfun(@(x) sum((b(d==x)).^2),w)./(4*w.^2); % total non-centrality for each eigenvalue

    b_sq_sum=accumarray(ic, b(nz).^2)';
    lambda=(b_sq_sum./(4*w.^2)); % total non-centrality for each eigenvalue
else
    w=d(nz);
    k=ones(size(w));
    lambda=b(nz).^2./(4*w.^2);
end

m=q0-dot(w,lambda);
s=norm(b(~nz));

% Optional 6th output: the eigen-structure of the standardized quadratic,
% shared by cdf_grad_norm_quad so it need not redo the eigendecomposition.
% These are the FULL (pre-merge) quantities: sqrt_v is the symmetric square
% root S=Sigma^{1/2}; V holds the eigenvectors and d the eigenvalues of
% S*Q2*S (its eigenvalues are the gx2 weights, with multiplicity the dof);
% b=V'*S*(2*Q2*mu+q1) as above.
if nargout>=6
    aux.S=sqrt_v;   % symmetric square root of the covariance, Sigma^{1/2}
    aux.V=R;        % eigenvectors of S*Q2*S (orthogonal)
    aux.d=d;        % row of eigenvalues of S*Q2*S (full, before merging)
    aux.b=b;        % row, V'*S*(2*Q2*mu + q1)
end