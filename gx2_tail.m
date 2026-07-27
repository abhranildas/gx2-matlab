function p=gx2_tail(x,w,k,l,s,m,varargin)

parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'x',@(x) isreal(x));
addRequired(parser,'w',@(x) isreal(x) && isrow(x));
addRequired(parser,'k',@(x) isreal(x) && isrow(x));
addRequired(parser,'l',@(x) isreal(x) && isrow(x));
addRequired(parser,'s',@(x) isreal(x) && isscalar(x));
addRequired(parser,'m',@(x) isreal(x) && isscalar(x));
addOptional(parser,'side','lower',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );
addParameter(parser,'output','cdf',@(x) strcmpi(x,'cdf') || strcmpi(x,'pdf') );
% addParameter(parser,'scale','linear',@(x) strcmpi(x,'linear') || strcmpi(x,'log') );

parse(parser,x,w,k,l,s,m,varargin{:});

% first merge into unique w's
[w,~,ic]=uniquetol(w); % unique non-zero eigenvalues
k=arrayfun(@(x)sum(k(ic==x)),1:numel(w)); % merged total dof's
l=arrayfun(@(x)sum(l(ic==x)),1:numel(w)); % merged total non-centralities

if strcmpi(parser.Results.side,'upper') % upper tail
    [w_max,max_idx]=max(w.*(w>0));
elseif strcmpi(parser.Results.side,'lower') % lower tail
    [w_max,max_idx]=min(w.*(w<0));
end

k_max=k(max_idx);
l_max=l(max_idx);

w_rest=w([1:max_idx-1, max_idx+1:end]);
k_rest=k([1:max_idx-1, max_idx+1:end]);
l_rest=l([1:max_idx-1, max_idx+1:end]);

a=exp(m/(2*w_max)+s^2/(8*w_max^2))*...
    prod(exp((l_rest.*w_rest)./(2*(w_max-w_rest)))./(1-w_rest/w_max).^(k_rest/2));

if strcmpi(parser.Results.output,'pdf')
    % if strcmpi(parser.Results.scale,'linear')
    p=a/abs(w_max)*ncx2pdf(x/w_max,k_max,l_max);
    % use log for tiny probabilities
    x_tiny=x(~p);

    % elseif strcmpi(parser.Results.scale,'log')
    if l_max % if non-central
        p_tiny=log10(a)-log10(abs(w_max))-log10(2*sqrt(2*pi))+(k_max-3)/4*log10(x_tiny/w_max)-(k_max-1)/4*log10(l_max)+(sqrt(l_max*x_tiny/w_max)-(l_max+x_tiny/w_max)/2)/log(10);
    else % if central
        p_tiny=log10(a)-log10(abs(w_max))-k_max/2*log10(2)-log10(gamma(k_max/2))+(k_max/2-1)*log10(x_tiny/w_max)-x_tiny/(2*w_max*log(10));
        % end
    end
elseif strcmpi(parser.Results.output,'cdf')
    % marcum=marcumq(sqrt(l_max),sqrt(x/w_max),k_max/2);
    p=a*ncx2cdf(x/w_max,k_max,l_max,'upper');

    % use log for tiny probabilities
    x_tiny=x(~p);
    if l_max % if non-central
        p_tiny=log10(a)-log10(l_max^((k_max-1)/4)*sqrt(2*pi))-(sqrt(x_tiny/w_max)-sqrt(l_max)).^2/(2*log(10))+(k_max-3)/4*log10(x_tiny/w_max);
    else % if central
        p_tiny=log10(a)+((k_max-2)/2)*log10(x_tiny/(2*w_max))-x_tiny/(2*w_max*log(10))-log10(gamma(k_max/2));
    end
end
p_tiny(p_tiny==-inf)=0; % set exactly 0 to 0
p(~p)=p_tiny;
if any(p<0)
    warning('Some output values are too small for double precision. Returning their log10 values, which are negative.')
end


% end