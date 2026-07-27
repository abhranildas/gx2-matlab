function [p,p_err]=gx2cdf_ray(x,w,k,l,s,m,varargin)

    parser=inputParser;
    parser.KeepUnmatched=true;
    addRequired(parser,'x',@(x) isreal(x));
    addRequired(parser,'w',@(x) isreal(x));
    addRequired(parser,'k',@(x) isreal(x));
    addRequired(parser,'l',@(x) isreal(x));
    addRequired(parser,'s',@(x) isreal(x) && isscalar(x));
    addRequired(parser,'m',@(x) isreal(x) && isscalar(x));
    addOptional(parser,'side','lower',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );

    parse(parser,x,w,k,l,s,m,varargin{:});
    side=parser.Results.side;

    y=x(:); % flatten input array x into a column vector of levels y

    % find standard quadratic form corresponding to the gx2:
    quad=gx2_to_norm_quad_params(w,k,l,s,m);

    dim=numel(quad.q1);
    mu=zeros(dim,1);
    v=eye(dim);

    % integrate with 'lower' for lower side.
    [p,p_err]=int_norm_ray(mu,v,quad,varargin{:},'side',side,'fun_level',y);

    % reshape flattened output arrays to shape of input x
    p=reshape(p,size(x));
    if ~isempty(p_err)
        p_err=reshape(p_err,size(x));
    end
end