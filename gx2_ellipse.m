function [p,p_rel_err]=gx2_ellipse(x,w,k,lambda,m,varargin)

parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'x',@(x) isreal(x));
addRequired(parser,'w',@(x) isreal(x) && isrow(x)  && (all(x>0)||all(x<0)) );
addRequired(parser,'k',@(x) isreal(x) && isrow(x));
addRequired(parser,'lambda',@(x) isreal(x) && isrow(x));
addRequired(parser,'m',@(x) isreal(x) && isscalar(x));
addOptional(parser,'side','lower',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );
addParameter(parser,'output','cdf',@(x) strcmpi(x,'cdf') || strcmpi(x,'pdf') );
addParameter(parser,'x_scale','linear',@(x) strcmpi(x,'linear') || strcmpi(x,'log') );

parse(parser,x,w,k,lambda,m,varargin{:});
side=parser.Results.side;
x_scale=parser.Results.x_scale;
w_pos=true;

% flatten x:
x_flat=x(:);

if all(w<0)
    w=-w; m=-m; w_pos=false;
    if strcmpi(x_scale,'linear')
        x_flat=-x_flat;
    end
end

% find ellipse parameters
ellipse_center=arrayfun(@(lambda,k) [sqrt(lambda) zeros(1,k-1)],lambda,k,'un',0);
ellipse_center=[ellipse_center{:}];
ellipse_weights=arrayfun(@(w,k) w*ones(1,k),w,k,'un',0);
ellipse_weights=[ellipse_weights{:}];

dim=sum(k);
% common factor to cdf and pdf (multinormal density at the ellipse center
% times the per-unit-x ellipse volume factor):
a=exp(-norm(ellipse_center)^2/2)/(2^(dim/2)*gamma(dim/2+1)*sqrt(prod(ellipse_weights)));

if strcmpi(x_scale,'linear')
    x_eff=max(x_flat-m,0); % effective x value from the tail
    if strcmpi(parser.Results.output,'cdf')
        p=a.*x_eff.^(dim/2);
        % flip if necessary
        if (w_pos && strcmpi(side,'upper')) || (~w_pos && strcmpi(side,'lower'))
            p=1-p;
        end
    elseif strcmpi(parser.Results.output,'pdf')
        p=(a*dim/2)*(x_flat-m).^(dim/2-1);
    end

else
    % compute log p for log x
    log10_x=x_flat;

    if strcmpi(parser.Results.output,'cdf')
        log10_p=dim/2*(log10_x-log10(2)) - norm(ellipse_center)^2/log(100) - log10(gamma(dim/2+1)) - sum(log10(ellipse_weights))/2;
    elseif strcmpi(parser.Results.output,'pdf')
        log10_p=(dim/2-1)*log10_x-(dim/2+1)*log10(2) +log10(dim) - norm(ellipse_center)^2/log(100) - log10(gamma(dim/2+1)) - sum(log10(ellipse_weights))/2;
    end
    p=log10_p;
end

% compute the lower and upper bounds for the estimate

if norm(ellipse_center) % if non-central
    if strcmpi(x_scale,'linear')
        radius_ratio=sqrt(x_eff)/sqrt(sum(ellipse_center.^2.*ellipse_weights));
        p_rel_err=norm(ellipse_center)^2*radius_ratio;
        % compute log (rel err) for log x
    else
        p_rel_err=2*log10(norm(ellipse_center))+log10_x/2-log10(sum(ellipse_center.^2.*ellipse_weights))/2;
    end
else % if central
    if strcmpi(x_scale,'linear')
        p_rel_err=x_eff/(2*min(ellipse_weights));
    else
        % compute log (rel err) for log x
        p_rel_err=log10_x-log10(2*min(ellipse_weights));
    end
end

% reshape outputs to input shape
p=reshape(p,size(x));
p_rel_err=reshape(p_rel_err,size(x));