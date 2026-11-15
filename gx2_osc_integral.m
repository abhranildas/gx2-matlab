function [I,failed]=gx2_osc_integral(H,c,AbsTol,RelTol)
% GX2_OSC_INTEGRAL Fourier-inversion integrals over v in (0,inf) of
% real(H(v).*exp(-1i*c*v)), shared by gx2_imhof and boundary_raw.
%
% H         amplitude: for a row vector v, returns a P-by-numel(v) complex
%           array (P=1 for gx2_imhof, one entry per derivative block for
%           boundary_raw). It is expressed in a dimensionless variable v in
%           which the integrand's mass lies at v~1.
% c         row vector of frequencies, one per point.
% I         P-by-numel(c) integrals.
% failed    1-by-numel(c), true where a point's integrals did not reach the
%           requested accuracy.
%
% The amplitude varies slowly, but can decay as slowly as a power of v (e.g.
% at s=0, and more so when some weights are tiny), so the integrand can keep
% oscillating far out, which a generic quadrature over (0,inf) cannot always
% handle. Several points are first integrated together in one shared plain
% quadrature over (0,inf), which is fast and handles most cases. If that
% does not converge, or there is only one point, each point is integrated
% separately. Within |c|<=1 (2 sd of m) the oscillation is slow over the
% integrand's mass, so a plain quadrature usually works and is more
% accurate, and is tried first. Otherwise, or if it fails, the integral is
% split (see fourier_split below) into a head integrated plainly and a tail
% integrated cycle by cycle, with the partial sums extrapolated by Wynn's
% epsilon algorithm, as QUADPACK's Fourier-integral routine QAWF does.

% Abhranil Das
% Center for Perceptual Systems, University of Texas at Austin
% Comments, questions, bugs to abhranil.das@utexas.edu

P=size(H(0),1);
% evaluation budget of the plain quadratures over (0,inf), about 200
% intervals as in the Python port (quad_vec's limit=200), so that a failed
% attempt is cheap (integral's own cap is 650 intervals, not settable)
MAX_PLAIN=3000;
nc=numel(c);
failed=false(1,nc);
if nc>1
    [I,converged]=quiet_integral(@(v) real(H(v).*exp(-1i*c*v)),0,inf,AbsTol,RelTol,MAX_PLAIN);
    if converged, return, end
end
I=nan(P,nc);
split=abs(c)>1;
for i=find(~split)
    % at c=0 there is no split to fall back on, so no budget either
    [I(:,i),converged]=quiet_integral(@(v) real(H(v).*exp(-1i*c(i)*v)),0,inf,AbsTol,RelTol,MAX_PLAIN/(c(i)~=0));
    failed(i)=~converged;
    split(i)=~converged && c(i)~=0;
end
% all points to split together; if the shared quadratures do not converge,
% redo them one at a time so that only the points that truly fail are
% flagged
if any(split)
    [I(:,split),failed(split),converged]=fourier_split(H,c(split),AbsTol,RelTol);
    if ~converged && nnz(split)>1
        for i=find(split)
            [I(:,i),failed(i)]=fourier_split(H,c(i),AbsTol,RelTol);
        end
    end
end
end

function [I,failed,converged]=fourier_split(H,c,AbsTol,RelTol)
% Integrals for a row vector of nonzero c; failed flags points whose
% integrals did not reach the requested accuracy, and converged is false if
% the shared quadratures did not.
%
% With c~=0 the integrand oscillates with period 2*pi/|c|. We integrate a
% head (0,v0) plainly, and the tail the way QUADPACK's QAWF does (MATLAB has
% no built-in equivalent): over successive cycles, each an odd number of
% half-periods long so that the cycle integrals alternate in sign, then
% extrapolate their partial sums with Wynn's epsilon algorithm. As in QAWF,
% a cycle is (2*floor(|c|)+1) half-periods, i.e. about 2*pi long for
% |c|>=1, and 50 cycles are used. The head is one period, but no longer than
% 2*pi, so that it covers the integrand's mass even for small c. The heads
% of all points, and all cycles of all points, are each integrated in one
% vectorized quadrature, over a variable t in (0,1) scaled to each point's
% interval.
n_cycles=50;
nc=numel(c);
v0=2*pi./max(abs(c),1);
[head,head_ok]=quiet_integral(@(t) real(H(v0*t).*exp(-1i*c.*v0*t)).*v0,0,1,AbsTol,RelTol);
P=size(head,1);
L=(2*floor(abs(c))+1)*pi./abs(c);
starts=v0+(0:n_cycles-1)'*L;      % n_cycles x nc
cycle=@(t) real(reshape(H(reshape(starts+t*L,1,[])),P,n_cycles,nc) ...
    .*reshape(exp(-1i*c.*(starts+t*L)),1,n_cycles,nc)).*reshape(L,1,1,nc);
[cycle_integrals,cycles_ok]=quiet_integral(cycle,0,1,AbsTol/n_cycles,RelTol);
converged=head_ok && cycles_ok;
I=nan(P,nc); failed=false(1,nc);
for i=1:nc
    for p=1:P
        [tail,err]=wynn_epsilon(cumsum(cycle_integrals(p,:,i)));
        I(p,i)=head(p,i)+tail;
        failed(i)=failed(i) || err>max(AbsTol,RelTol*abs(I(p,i)));
    end
    failed(i)=failed(i) || ~converged;
end
end

function [val,err]=wynn_epsilon(S)
% Wynn's epsilon algorithm: the extrapolated limit of the sequence of
% partial sums S, and an error estimate. Each even column of the epsilon
% table is a sequence of estimates of the limit; we return the last entry of
% the even column that changed least from the previous even column.
n=numel(S);
val=S(end);
err=abs(S(end)-S(max(end-1,1)));
prev=zeros(1,n+1);         % column -1
cur=S(:)';                 % column 0
last_even=cur(end);
for col=1:n-1
    nxt=prev(2:numel(cur))+1./diff(cur);
    if ~all(isfinite(nxt)), break, end % sequence already converged exactly
    prev=cur; cur=nxt;
    if mod(col,2)==0
        e=abs(cur(end)-last_even);
        if e<err, val=cur(end); err=e; end
        last_even=cur(end);
    end
end
end

function [q,converged]=quiet_integral(f,a,b,AbsTol,RelTol,max_evals)
% integral (array-valued) with its convergence warnings captured instead of
% shown: converged is false if integral warned that it did not reach the
% requested accuracy, or if it was stopped after max_evals evaluations of f
% (optional, default no limit; q is then nan).
if nargin<6, max_evals=inf; end
n_evals=0;
ids={'MATLAB:integral:MaxIntervalCountReached','MATLAB:integral:MinStepSize',...
    'MATLAB:integral:NonFiniteValue'};
state=warning;
cellfun(@(id) warning('off',id),ids);
[last_msg,last_id]=lastwarn;
lastwarn('');
try
    q=integral(@counted,a,b,'AbsTol',AbsTol,'RelTol',RelTol,'ArrayValued',true);
    [~,id]=lastwarn;
    converged=~any(strcmp(id,ids));
catch err
    if ~strcmp(err.identifier,'gx2:evalBudget'), warning(state); rethrow(err); end
    q=nan(size(f(a))); converged=false; id='';
end
warning(state);
if isempty(id), lastwarn(last_msg,last_id); end

    function y=counted(v)
        n_evals=n_evals+1;
        if n_evals>max_evals, error('gx2:evalBudget','evaluation budget exhausted'); end
        y=f(v);
    end
end
