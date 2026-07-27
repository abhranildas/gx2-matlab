function [p,p_err]=gx2_ruben_eval(coeffs,x,m,varargin)

% GX2_RUBEN_EVAL Evaluate Ruben's series at x (offset m) from coefficients
% already computed by gx2_ruben_coeffs.
%
% Abhranil Das
% Center for Perceptual Systems, University of Texas at Austin
% Comments, questions, bugs to abhranil.das@utexas.edu
%
% Usage:
% p=gx2_ruben_eval(coeffs,x,m)
% p=gx2_ruben_eval(coeffs,x,m,'upper','output','pdf','nx',1)
%
% Required inputs:
% coeffs    struct returned by gx2_ruben_coeffs
% x         array of points to evaluate at
% m         scalar offset
%
% Optional positional/name-value inputs:
% side      'lower' (default) or 'upper'
% output    'cdf' (default) or 'pdf'
% nx        x-derivative order (pdf output only). Default=0.
%
% Outputs:
% p         cdf or pdf values, same size as x
% p_err     truncation-error indicator (leftover series mass times the next
%           central-chi-square factor)
%
% See also:
% gx2_ruben_coeffs, gx2_ruben

parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'coeffs',@isstruct);
addRequired(parser,'x',@(x) isreal(x));
addRequired(parser,'m',@(x) isreal(x) && isscalar(x));
addOptional(parser,'side','lower',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );
addParameter(parser,'output','cdf',@(x) strcmpi(x,'cdf') || strcmpi(x,'pdf') );
addParameter(parser,'nx',0,@(x) isscalar(x) && (x>=0) && (x==round(x)));

parse(parser,coeffs,x,m,varargin{:});
side=parser.Results.side;
nx=parser.Results.nx;
if nx>0 && ~strcmpi(parser.Results.output,'pdf')
    error('The x-derivative order ''nx'' is only defined for the ''pdf'' output.')
end

a=coeffs.a; N=coeffs.N; beta=coeffs.beta; M=coeffs.M; w_pos=coeffs.w_pos;

% flatten x, and flip into the coefficients' all-positive-weight frame
x_flat=x(:);
if ~w_pos
    x_flat=-x_flat; m=-m;
end

% compute the central chi-squared integrals (only the terms actually used)
[x_grid,k_grid]=meshgrid((x_flat-m)/beta,M:2:M+2*(N-1));
if strcmpi(parser.Results.output,'cdf')
    if (w_pos && strcmpi(side,'upper')) || (~w_pos && strcmpi(side,'lower'))
        % upper tail
        F=chi2cdf(x_grid,k_grid,'upper');
    else
        F=chi2cdf(x_grid,k_grid);
    end
elseif strcmpi(parser.Results.output,'pdf')
    F=chi2pdf_nderiv(x_grid,k_grid,nx);   % nx-th y-derivative of the chi2 density
end

% compute the integral
p=a'*F;

if strcmpi(parser.Results.output,'pdf')
    % each x-derivative brings a factor 1/beta from y=(x-m)/beta; the flipped
    % (all-negative-weight) frame contributes a factor (-1)^nx.
    p=p/beta^(nx+1);
    if ~w_pos, p=p*(-1)^nx; end
end

% truncation-error indicator: the leftover series mass (now negligible unless
% the n_ruben cap was hit) times the next central-chi-square factor
p_err=(1-sum(a))*chi2cdf((x_flat-m)/beta,M+2*N);

% reshape outputs to input shape
p=reshape(p,size(x));
p_err=reshape(p_err,size(x));

end

% ---------------------------------------------------------------------------
function gd=chi2pdf_nderiv(y,nu,n)
% n-th derivative in y of the central chi-square density g_nu(y). Uses the
% closed form  g_nu^(n)(y) = g_nu(y) * sum_{j=0}^n C(n,j)(-1/2)^{n-j} (a)_j y^{-j},
% where a=nu/2-1 and (a)_j = a(a-1)...(a-j+1) is the falling factorial. This
% is exact for any nu>0 (no negative-dof chi-square ever appears). The
% derivative vanishes on the support edge y<=0.
gd=chi2pdf(y,nu);
if n==0, return; end
a=nu/2-1;
poly=zeros(size(gd));
for j=0:n
    ff=ones(size(a));
    for l=0:j-1
        ff=ff.*(a-l);       % falling factorial (a)_j, elementwise
    end
    poly=poly+nchoosek(n,j)*(-0.5)^(n-j).*ff.*(y.^(-j));
end
gd=gd.*poly;
gd(y<=0)=0;
end
