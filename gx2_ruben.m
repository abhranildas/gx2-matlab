function [p,p_err]=gx2_ruben(x,w,k,lambda,m,varargin)

parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'x',@(x) isreal(x));
addRequired(parser,'w',@(x) isreal(x) && isrow(x)  && (all(x>0)||all(x<0)) );
addRequired(parser,'k',@(x) isreal(x) && isrow(x));
addRequired(parser,'lambda',@(x) isreal(x) && isrow(x));
addRequired(parser,'m',@(x) isreal(x) && isscalar(x));
addOptional(parser,'side','lower',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );
addParameter(parser,'output','cdf',@(x) strcmpi(x,'cdf') || strcmpi(x,'pdf') );
addParameter(parser,'n_ruben',1e3,@(x) isscalar(x) && (x>0) && (x==round(x)));
addParameter(parser,'nx',0,@(x) isscalar(x) && (x>=0) && (x==round(x)));

parse(parser,x,w,k,lambda,m,varargin{:});
side=parser.Results.side;
n_ruben=parser.Results.n_ruben;
nx=parser.Results.nx;   % x-derivative order (pdf output only): 0 gives the pdf
if nx>0 && ~strcmpi(parser.Results.output,'pdf')
    error('The x-derivative order ''nx'' is only defined for the ''pdf'' output.')
end

% flatten x:
x_flat=x(:);

w_pos=true;
if all(w<0)
    w=-w; x_flat=-x_flat; m=-m; w_pos=false;
end
beta=0.90625*min(w);
M=sum(k);
n=(1:n_ruben-1)';

% compute the g's
g=sum(k.*(1-beta./w).^n,2)+ beta*n.*((1-beta./w).^(n-1))*(lambda./w)';

% compute the expansion coefficients, stopping once the leftover series mass
% is negligible. The a_j are nonnegative and sum to 1, so the tail mass
% 1-sum(a_{1:N}) both bounds the truncation error and decreases monotonically.
% The stop uses only this cheap coefficient recurrence -- not the chi-square
% grid below -- so the term count N is fixed in a single pass, and the
% expensive evaluation is then done once at that reduced size. n_ruben is the
% safety cap; most cases converge in ~10^2 terms well under it.
masstol=1e-14;
a=nan(n_ruben,1);
a(1)=sqrt(exp(-sum(lambda))*beta^M*prod(w.^(-k)));
if a(1)<realmin
    error('Underflow error: some series coefficients are smaller than machine precision.')
end
cum=a(1); N=n_ruben;
for j=1:n_ruben-1
    a(j+1)=dot(flip(g(1:j)),a(1:j))/(2*j);
    cum=cum+a(j+1);
    if 1-cum<masstol, N=j+1; break; end
end
a=a(1:N);

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

if strcmpi(parser.Results.output,'cdf')
    % flip if necessary
    if (w_pos && strcmpi(side,'upper')) || (~w_pos && strcmpi(side,'lower'))
        % abar=1-sum(a);
        % p=p+abar;
    end
elseif strcmpi(parser.Results.output,'pdf')
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