function [p,errflag]=gx2_imhof(x,w,k,l,s,m,varargin)
% Imhof-Davies method for the cdf/pdf of a generalized chi-square, and for
% the exact (non-finite-differenced) parameter-derivative integrands used by
% cdf_grad_gx2 -- see gx2_imhof_integrand for the 'output' modes.
%
% The result is an integral over u in (0,inf) of an integrand that
% oscillates with period 4*pi/|x-m| under an amplitude that can decay as
% slowly as a power of u (e.g. when some weights are tiny), which a generic
% quadrature over (0,inf) cannot always handle. With 'precision' 'basic',
% several x points are first integrated together in one shared plain
% quadrature over (0,inf), which is fast and handles most cases. If that does
% not converge, or there is only one point, each point is integrated
% separately: within 2 sd of m with a plain quadrature first, and otherwise,
% or if that fails, split into a
% head over the first period (capped at 2*pi/sd in u), integrated plainly,
% and a tail integrated cycle by cycle with the partial sums extrapolated by
% Wynn's epsilon algorithm, as QUADPACK's Fourier-integral routine QAWF does
% (see gx2_osc_integral); this is vectorized over all x points. All
% quadratures run in the dimensionless variable sd*u, so results do not
% depend on the units of x.
%
% AbsTol and RelTol are tolerances on the integral itself, not on the
% returned probability. For the cdf the integral is of order pi/2, so the
% defaults allow an absolute error of up to about 5e-7 in p (typically far
% less); for small tail probabilities, tighten them (e.g. AbsTol=1e-14,
% RelTol=1e-12) for tail-relative accuracy. Probabilities below about 1e-14
% are lost to cancellation in 0.5-/+integral/pi at any tolerance; use the
% 'tail' or 'ray' methods there.
%
% 'precision' 'vpa' uses a plain variable-precision quadrature over (0,inf),
% without the treatment of the oscillatory tail above, so it can be
% inaccurate (and slow) when the integrand decays slowly far from m.
%
% errflag is true where p is unreliable: either the quadrature did not reach
% the requested accuracy (warning gx2:imhofNoConvergence), or a probability
% fell outside [0,1] and was clipped (warning gx2:imhofClip). It is not set
% by the 'vpa' path's quadrature, which reports no convergence status.

parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'x',@(x) isreal(x));
addRequired(parser,'w',@(x) isreal(x) && isrow(x));
addRequired(parser,'k',@(x) isreal(x) && isrow(x));
addRequired(parser,'l',@(x) isreal(x) && isrow(x));
addRequired(parser,'s',@(x) isreal(x) && isscalar(x));
addRequired(parser,'m',@(x) isreal(x) && isscalar(x));
addOptional(parser,'side','lower',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );
addParameter(parser,'output','cdf',@(x) any(strcmpi(x,{'cdf','pdf','dens','k_deriv','kk_deriv'})) );
addParameter(parser,'idx',[],@(x) isempty(x) || (isvector(x) && all(x==round(x)) && all(x>=1)));
addParameter(parser,'nx',0,@(x) isscalar(x) && x==round(x) && x>=0);
addParameter(parser,'precision','basic',@(x) strcmpi(x,'basic')||strcmpi(x,'vpa'));
addParameter(parser,'AbsTol',1e-10,@(x) isreal(x) && isscalar(x) && (x>=0));
addParameter(parser,'RelTol',1e-6,@(x) isreal(x) && isscalar(x) && (x>=0));

parse(parser,x,w,k,l,s,m,varargin{:});
output=parser.Results.output;
idx=parser.Results.idx;
nx=parser.Results.nx;
side=parser.Results.side;
AbsTol=parser.Results.AbsTol;
RelTol=parser.Results.RelTol;

% Integrate in the dimensionless variable v=sd*u, in which the integrand's
% mass lies at v~1 and it oscillates at rate (x-m)/(2*sd). The quadrature
% routines have built-in unit scales (the map of (0,inf), the cycle length),
% so without this their accuracy would depend on the units of x.
[~,v]=gx2stat(w,k,l,s,m);
sd=sqrt(v);
if sd==0, sd=1; end
cs=(x(:)'-m)/(2*sd);
failed=false(size(cs));

% compute the integral
if strcmpi(parser.Results.precision,'basic')
    H=@(v) gx2_imhof_amplitude(v/sd,w',k',l',s,output,idx,nx)/sd;
    [imhof_integral,failed]=gx2_osc_integral(H,cs,AbsTol,RelTol);
    if any(failed)
        warning('gx2:imhofNoConvergence',...
            'Imhof integral did not reach the requested accuracy at %d of %d point(s); these are flagged in errflag.',...
            nnz(failed),numel(failed))
    end
elseif strcmpi(parser.Results.precision,'vpa')
    % plain variable-precision quadrature over (0,inf). Unlike the 'basic'
    % path (and the Python port's 'vpa', which uses mpmath.quadosc), this
    % has no special treatment of the oscillatory tail, as MATLAB has no
    % variable-precision oscillatory quadrature, so it can be inaccurate when
    % the integrand decays slowly (e.g. tiny weights) far from m.
    syms u
    imhof_integral=arrayfun(@(x) vpaintegral(@(u) gx2_imhof_integrand(u,x,w',k',l',s,m,output,idx,nx),...
        u,0,inf,'AbsTol',AbsTol,'RelTol',RelTol),x(:)');
end
imhof_integral=reshape(imhof_integral,size(x));
failed=reshape(failed,size(x));

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
else
    % signed derivative outputs (no probability clipping):
    %   'dens'      x-derivatives of the pdf, normalized by 1/(2*pi)
    %   'k_deriv'   d/dk of the cdf (and its x-derivatives), normalized by 1/pi
    %   'kk_deriv'  d^2/(dk dk) of the cdf, normalized by 1/pi
    if strcmpi(output,'dens')
        p=imhof_integral/(2*pi);
    else
        p=imhof_integral/pi;
    end
    if isa(p,'sym')
        p=double(vpa(p));
    end
    errflag = false(size(p));
end

if any(errflag(:))
    warning('gx2:imhofClip','Imhof method output(s) too close to limit to compute exactly, so clipping. Check the flag output, and try stricter tolerances.')
    p=max(p,0);
end
errflag=errflag|failed;

end
