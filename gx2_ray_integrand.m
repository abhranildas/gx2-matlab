function [p_rays,p_tiny_sum,sym_idx]=gx2_ray_integrand(x,n_z,quad,varargin)
% return the differential probability or probability density on each ray
% that is integrated across rays

% parse inputs
parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'x',@isnumeric); % points at which to find the cdf/pdf
addRequired(parser,'n_z',@isnumeric);
addRequired(parser,'quad',@isstruct);
addOptional(parser,'side','upper',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );
addParameter(parser,'output','prob'); % probability or probability density
addParameter(parser,'precision','log',@(x) strcmpi(x,'basic')||strcmpi(x,'log')||strcmpi(x,'vpa'));

parse(parser,x,n_z,quad,varargin{:});

output=parser.Results.output;
side=parser.Results.side;
precision=parser.Results.precision;

dim=numel(quad.q1);
n_levels=numel(x);

% find the quadratic coefficients across all rays
q2=dot(n_z,quad.q2*n_z,1);
q1=quad.q1'*n_z;
q0=quad.q0;

% discriminant of the quadratic across rays and levels
delta2=q1.^2-4*q2.*(q0-x); % delta^2
root_exists=delta2>0; % levels where linear or quadratic roots exist
quad_root_exists=root_exists & q2; % levels where quadratic roots exist
delta=nan(size(delta2));
delta(quad_root_exists)=sqrt(delta2(quad_root_exists)); % populate only with quad delta for now

% linear_root_exists=repmat(~q2 & q1, [numel(x) 1]); % levels where linear roots exist

% sorted roots across rays
z=(-q1+cat(3,-1,1).*delta)./(2*abs(q2)); % quadratic roots where q2 ~= 0
if nnz(~q2)
    z(:,~q2,1)=-(q0-x)./q1(~q2); % linear roots where q2=0
end

if strcmpi(output,'prob')
    init_sign_rays=sign(4*sign(q2)-2*sign(q1)+sign(q0-x));
    [Phibar_big,Phibar_small]=Phibar_ray_split(z,dim);
    p_rays_big=init_sign_rays+1+init_sign_rays.*(Phibar_big(:,:,2)-Phibar_big(:,:,1));
    p_rays_small=init_sign_rays.*(Phibar_small(:,:,2)-Phibar_small(:,:,1));
    if strcmpi(side,'upper')
        p_rays=p_rays_big+p_rays_small;
    elseif strcmpi(side,'lower')
        p_rays=2-p_rays_big-p_rays_small;
    end
elseif strcmpi(output,'prob_dens')
    % phi_ray at each root
    sum_phi=sum(phi_ray(z,dim),3,'omitnan');

    % quadratic slope at each root
    quad_slope=nan(size(delta2));
    quad_slope(root_exists)=sqrt(delta2(root_exists)); % slope is the same formula at quad and linear roots

    % divide phi_ray by quad slope
    p_rays=sum_phi./quad_slope;

    p_rays(isnan(p_rays))=0;
end

% now compute tiny probabilities or densities < realmin

tiny_probs=root_exists&(~p_rays); % cases where root exists but computed prob. or density is 0

if strcmpi(precision,'basic')
    if nnz(tiny_probs)
        if strcmpi(output,'prob')
            warning("%.1f%% of rays contain probabilities less than realmin=1e-308, returning 0. Set 'precision' to 'log' or 'vpa' to compute these.",100*mean(tiny_probs,'all'))
        elseif strcmpi(output,'prob_dens')
            warning("%.1f%% of rays contain probability densities less than realmin=1e-308, returning 0. Set 'precision' to 'log' or 'vpa' to compute these.",100*mean(tiny_probs,'all'))
        end
    end
elseif strcmpi(precision,'log')
    % find the roots where value was tiny.
    % for prob, these are where Phibar_small was tiny, or where they weren't tiny
    % individually but cancelled across the two roots
    if nnz(tiny_probs) % only do the log-domain work if any tiny values exist
    tiny_probs=repmat(tiny_probs,[1 1 2]);
    z_tiny=nan(size(z));
    z_tiny(tiny_probs)=z(tiny_probs);

    if strcmpi(output,'prob')
        % sign of contribution from each root
        z_tiny_signs=init_sign_rays.*sign(z_tiny).*cat(3,-1,1);
        if strcmpi(side,'lower')
            z_tiny_signs=-z_tiny_signs;
        end

        % log of tiny Phibar smalls:
        p_tiny=nan(size(z_tiny));

        % lower tail prob. of the chi distribution < realmin,
        % z<sqrt(chi2inv(realmin,dim)), but we'll just use the median as
        % an easy cutoff
        z_med=sqrt(chi2inv(.5,dim));
        z_tiny_lo=abs(z_tiny)<z_med;
        p_tiny(z_tiny_lo)=dim*log10(abs(z_tiny(z_tiny_lo)))-log10(gamma(dim/2)*2^(dim/2-1)*dim);

        % upper tail prob. of the chi distribution < realmin:
        % z_hi_thresh=sqrt(fzero(@(x) chi2cdf(x,dim,'upper') - realmin, [0, 1e6]));
        z_tiny_hi=abs(z_tiny)>z_med;
        p_tiny(z_tiny_hi)=(dim-2)*log10(abs(z_tiny(z_tiny_hi)))-z_tiny(z_tiny_hi).^2/(2*log(10))-log10(gamma(dim/2)*2^(dim/2-1));

        % combine magnitude and sign
        signed_p_tiny=p_tiny.*z_tiny_signs;

        % now do signed log sum exp
        p_tiny_sum=signed_log_sum_exp(signed_p_tiny,[2 3]);

    elseif strcmpi(output,'prob_dens')
        p_tiny=(dim-1)*log10(abs(z_tiny))-z_tiny.^2/(2*log(10))-log10(gamma(dim/2)*2^(dim/2-1));
        p_tiny(isnan(p_tiny))=-inf;
        % now divide by quad slope and sum using log sum exp
        p_tiny_sum=signed_log_sum_exp(p_tiny,3)-log10(2*quad_slope);
        p_tiny_sum=signed_log_sum_exp(p_tiny_sum,2);
    end
    else % no tiny values exist; nothing to add in the log domain
        p_tiny_sum=-inf(n_levels,1);
    end



elseif strcmpi(precision,'vpa')
    p_tiny_sum=num2cell(zeros(n_levels,1));
    sym_idx=[];
    if nnz(tiny_probs)
        sym_idx=any(tiny_probs,2); % levels at which at least one ray is sym.
        % group roots into a cell array so we can use cellfun:
        z_cell=permute(num2cell(permute(z,[1 3 2]),2),[1 3 2]);
        % remove nan roots:
        z_cell=cellfun(@(z_each) z_each(~isnan(z_each)),z_cell,'un',0);
        if strcmpi(output,'prob')
            p_tiny_sum(sym_idx)=sym2cell(sum(cellfun(@(init_sign_ray,z_ray) prob_ray_sym(init_sign_ray,z_ray,dim,side),num2cell(init_sign_rays(sym_idx,:)),z_cell(sym_idx,:)),2));
        elseif strcmpi(output,'prob_dens')
            p_sym=cellfun(@(z_ray,quad_slope_ray) sum(phi_ray(sym(z_ray),dim))/quad_slope_ray,z_cell(sym_idx,:),num2cell(quad_slope(sym_idx,:)));
            p_tiny_sum(sym_idx)=sym2cell(sum(p_sym,2,'omitnan'));
        end
    end
end