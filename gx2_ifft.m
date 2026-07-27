function [p,xgrid]=gx2_ifft(x,w,k,l,s,m,varargin)

parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'x',@(x) isreal(x) || strcmpi(x,'full'));
addRequired(parser,'w',@(x) isreal(x) && isrow(x));
addRequired(parser,'k',@(x) isreal(x) && isrow(x));
addRequired(parser,'l',@(x) isreal(x) && isrow(x));
addRequired(parser,'s',@(x) isreal(x) && isscalar(x));
addRequired(parser,'m',@(x) isreal(x) && isscalar(x));
addOptional(parser,'side','lower',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );
addParameter(parser,'output','cdf',@(x) strcmpi(x,'cdf') || strcmpi(x,'pdf') );
% range of x over which to compute ifft:
if strcmpi(x,'full')
    % default span is upto 100 sd from mean
    [mu,v]=gx2stat(w,k,l,s,m);
    span=max(abs(mu+[-1 1]*100*sqrt(v)));
else
    % default span is 100* input span
    % span=100*max(abs(x(:)));
    span=1e5;
end
addParameter(parser,'span',span,@(x) isscalar(x) && (x>0));
% number of grid points for ifft
addParameter(parser,'n_grid',1e6+1,@(x) isscalar(x) && (x>0));
addParameter(parser,'ft_type','cft');
parse(parser,x,w,k,l,s,m,varargin{:});
span=parser.Results.span;
n_grid=round(parser.Results.n_grid);

% make n_grid odd
if ~mod(n_grid,2)
    n_grid=n_grid+1;
end

n=(n_grid-1)/2;
idx=-n:n;
dx=span/n;

if strcmpi(x,'full')
    x_mid=0;
else
    x_mid=(min(x(:))+max(x(:)))/2;
end
xgrid=x_mid+idx*dx; % grid of x over which this cdf/pdf is computed

if strcmpi(parser.Results.ft_type,'dft')
    % compute non-central and normal pdf's over this span
    w_idx=find(w); % non-zero w indices
    ncpdfs=nan(nnz(w),length(xgrid));
    for i=1:nnz(w)
        if i==1 % add the offset m to the first term
            pdf=gx2pdf(xgrid,w(w_idx(i)),k(w_idx(i)),l(w_idx(i)),0,m);
        else
            pdf=gx2pdf(xgrid,w(w_idx(i)),k(w_idx(i)),l(w_idx(i)),0,0);
        end
        pdf(isinf(pdf))=max(pdf(~isinf(pdf)));
        ncpdfs(i,:)=pdf;
    end
    if s
        ncpdfs=[ncpdfs;normpdf(xgrid,0,abs(s))];
    end
    phi=prod(fft(ifftshift(ncpdfs,2),[],2),1);
    p=fftshift(ifft(phi));

    p=p/(sum(p)*dx); % normalize
elseif strcmpi(parser.Results.ft_type,'cft')
    dt=2*pi/(n_grid*dx);
    t=idx*dt; % grid of points to evaluate characteristic function
    phi=gx2char(-t,w,k,l,s,m);

    if strcmpi(parser.Results.output,'pdf')
        phi=phi.*exp(1i*x_mid*dt*idx);
        p=fftshift(ifft(ifftshift(phi)))/dx;
    elseif strcmpi(parser.Results.output,'cdf')
        phi=phi./(1i*t).*exp(1i*x_mid*dt*idx);
        phi(isinf(phi))=0;
        p=0.5+fftshift(ifft(ifftshift(phi)))/dx;
    end
end

if strcmpi(parser.Results.output,'cdf') && strcmpi(parser.Results.side,'upper')
    p=1-p;
end


if ~strcmpi(x,'full')
    F=griddedInterpolant(xgrid,p);
    p=F(x);
end

p=max(p,0);