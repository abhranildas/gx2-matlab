function logp=log_gx2cdf(x,w,k,l,s,m,varargin)
p=gx2cdf(x,w,k,l,s,m,varargin{:});

if p<=0
    logp=p;
elseif p>0
    logp=log10(p);
end
