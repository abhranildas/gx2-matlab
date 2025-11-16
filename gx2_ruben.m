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
addParameter(parser,'n_ruben',1e3,@(x) ismember(x,1:x));

parse(parser,x,w,k,lambda,m,varargin{:});
side=parser.Results.side;
n_ruben=parser.Results.n_ruben;

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

% compute the expansion coefficients
a=nan(n_ruben,1);
a(1)=sqrt(exp(-sum(lambda))*beta^M*prod(w.^(-k)));
if a(1)<realmin
    error('Underflow error: some series coefficients are smaller than machine precision.')
end
for j=1:n_ruben-1
    a(j+1)=dot(flip(g(1:j)),a(1:j))/(2*j);
end

% a=a/sum(a);

% compute the central chi-squared integrals
[x_grid,k_grid]=meshgrid((x_flat-m)/beta,M:2:M+2*(n_ruben-1));
if strcmpi(parser.Results.output,'cdf')
    if (w_pos && strcmpi(side,'upper')) || (~w_pos && strcmpi(side,'lower'))
        % upper tail
        F=arrayfun(@(x,k) chi2cdf(x,k,'upper'),x_grid,k_grid);
    else
        F=arrayfun(@(x,k) chi2cdf(x,k),x_grid,k_grid);
    end
elseif strcmpi(parser.Results.output,'pdf')
    F=arrayfun(@(x,k) chi2pdf(x,k),x_grid,k_grid);
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
    p=p/beta;
end

% compute the truncation error
p_err=(1-sum(a))*chi2cdf((x_flat-m)/beta,M+2*n_ruben);

% reshape outputs to input shape
p=reshape(p,size(x));
p_err=reshape(p_err,size(x));