function [p,errflag]=gx2_imhof(x,w,k,lambda,s,m,varargin)

parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'x',@(x) isreal(x));
addRequired(parser,'w',@(x) isreal(x) && isrow(x));
addRequired(parser,'k',@(x) isreal(x) && isrow(x));
addRequired(parser,'lambda',@(x) isreal(x) && isrow(x));
addRequired(parser,'s',@(x) isreal(x) && isscalar(x));
addRequired(parser,'m',@(x) isreal(x) && isscalar(x));
addOptional(parser,'side','lower',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );
addParameter(parser,'output','cdf',@(x) strcmpi(x,'cdf') || strcmpi(x,'pdf') );
addParameter(parser,'precision','basic',@(x) strcmpi(x,'basic')||strcmpi(x,'vpa'));
addParameter(parser,'AbsTol',1e-10,@(x) isreal(x) && isscalar(x) && (x>=0));
addParameter(parser,'RelTol',1e-6,@(x) isreal(x) && isscalar(x) && (x>=0));

parse(parser,x,w,k,lambda,s,m,varargin{:});
output=parser.Results.output;
side=parser.Results.side;
AbsTol=parser.Results.AbsTol;
RelTol=parser.Results.RelTol;

% compute the integral
if strcmpi(parser.Results.precision,'basic')
    % integrate over all x points in one adaptive quadrature: for each
    % quadrature node u, the x-independent parts of the integrand (theta's
    % sum over terms, and rho) are computed once and shared across all x.
    imhof_integral=integral(@(u) gx2_imhof_integrand(u,x(:)',w',k',lambda',s,m,output),...
        0,inf,'AbsTol',AbsTol,'RelTol',RelTol,'ArrayValued',true);
    imhof_integral=reshape(imhof_integral,size(x));
elseif strcmpi(parser.Results.precision,'vpa')
    syms u
    imhof_integral=arrayfun(@(x) vpaintegral(@(u) gx2_imhof_integrand(u,x,w',k',lambda',s,m,output),...
        u,0,inf,'AbsTol',AbsTol,'RelTol',RelTol),x);
end

if strcmpi(output,'cdf')
    if strcmpi(side,'lower')
        p=0.5-imhof_integral/pi;
    elseif strcmpi(side,'upper')
        p=0.5+imhof_integral/pi;
    end
    if isa(p,'sym')
        p=double(vpa(p));
    end
    errflag = p<0 | p>1;
    p=min(p,1);
elseif strcmpi(output,'pdf')
    p=imhof_integral/(2*pi);
    if isa(p,'sym')
        p=double(vpa(p));
    end
    errflag = p<0;
end

if any(errflag)
    warning('Imhof method output(s) too close to limit to compute exactly, so clipping. Check the flag output, and try stricter tolerances.')
    p=max(p,0);
end

end