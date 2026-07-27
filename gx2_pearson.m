function p=gx2_pearson(x,w,k,l,s,m,varargin)

parser = inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'x',@(x) isreal(x));
addRequired(parser,'w',@(x) isreal(x) && isrow(x));
addRequired(parser,'k',@(x) isreal(x) && isrow(x));
addRequired(parser,'l',@(x) isreal(x) && isrow(x));
addRequired(parser,'s',@(x) isreal(x) && isrow(x));
addRequired(parser,'m',@(x) isreal(x) && isscalar(x));
addOptional(parser,'side','lower',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );
addParameter(parser,'output','cdf',@(x) strcmpi(x,'cdf') || strcmpi(x,'pdf') );

parse(parser,x,w,k,l,s,m,varargin{:});
side=parser.Results.side;

% first 3 central moments:
mu(1)=sum(w.*(k+l))+m;
mu(2)=2*sum(w.^2.*(k+2*l))+s^2;
mu(3)=8*sum(w.^3.*(k+3*l));

h=8*mu(2)^3/mu(3)^2;

if mu(3)>0
    y=(x-mu(1))*sqrt(2*h/mu(2))+h;
    if strcmpi(parser.Results.output,'cdf')
        if strcmpi(side,'lower')
            p=chi2cdf(y,h);
        elseif strcmpi(side,'upper')
            p=chi2cdf(y,h,'upper');
        end
    elseif strcmpi(parser.Results.output,'pdf')
        p=sqrt(2*h/mu(2))*chi2pdf(y,h);
    end
else
    mu(1)=-mu(1);
    mu(3)=-mu(3);
    x=-x;
    y=(x-mu(1))*sqrt(2*h/mu(2))+h;
    if strcmpi(parser.Results.output,'cdf')
        if strcmpi(side,'lower')
            p=chi2cdf(y,h,'upper');
        elseif strcmpi(side,'upper')
            p=chi2cdf(y,h);
        end
    elseif strcmpi(parser.Results.output,'pdf')
        p=sqrt(2*h/mu(2))*chi2pdf(y,h);
    end
end

end