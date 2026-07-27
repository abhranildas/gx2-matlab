function [f,f_err]=gx2pdf_ray(x,w,k,l,s,m,varargin)

parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'x',@(x) isreal(x));
addRequired(parser,'w',@(x) isreal(x) && isrow(x));
addRequired(parser,'k',@(x) isreal(x) && isrow(x));
addRequired(parser,'l',@(x) isreal(x) && isrow(x));
addRequired(parser,'s',@(x) isreal(x) && isscalar(x));
addRequired(parser,'m',@(x) isreal(x) && isscalar(x));
addParameter(parser,'n_rays',1e3);

parse(parser,x,w,k,l,s,m,varargin{:});

y=x(:); % flatten input array x into a column vector of levels y

% find standard quadratic form corresponding to the gx2:
quad=gx2_to_norm_quad_params(w,k,l,s,m);
dim=numel(quad.q1);
mu=zeros(dim,1);
v=eye(dim);

% call int_norm_ray to integrate
[f,f_err]=int_norm_ray(mu,v,quad,'fun_level',y,'output','prob_dens',varargin{:});

% reshape flattened array to shape of input x
f=reshape(f,size(x));
if ~isempty(f_err)
    f_err=reshape(f_err,size(x));
end

end