function bench_cdf_grad_gx2(cases)
% BENCH_CDF_GRAD_GX2 Benchmark cdf_grad_gx2 against central finite
% differences of gx2cdf, for Table 2.1.2 of derivatives/gx2_derivatives.md.
% Same protocol as gx2-py/tests/bench_cdf_grad_gx2.py: ground truth is
% cdf_grad_gx2 at AbsTol=1e-13, RelTol=1e-11; FD is central differences of
% gx2cdf at default tolerances, swept over the step h, keeping the best
% error; errors are max-abs over entries, with and without the k entries.
% Usage: bench_cdf_grad_gx2        % cases 1:3
%        bench_cdf_grad_gx2(3)

if nargin<1, cases=1:3; end
C={struct('w',[1 -2],'k',[2 3],'l',[1 .5],'s',1.5,'m',1,'x',2), ...
   struct('w',[1 -2 3],'k',[1 2 1],'l',[.5 1 .5],'s',1,'m',2,'x',5), ...
   struct('w',[1 -1.5],'k',[1 1],'l',[.5 .3],'s',0,'m',0,'x',.5)};
tight={'AbsTol',1e-13,'RelTol',1e-11};
steps=logspace(-6,-1,11);

for cs=cases
    c=C{cs}; n=numel(c.w); x=c.x;
    p=[c.w c.k c.l c.s c.m]; P=numel(p);
    noK=true(1,P); noK(n+1:2*n)=false; noK2=noK'&noK;
    F=@(q) gx2cdf(x,q(1:n),q(n+1:2*n),q(2*n+1:3*n),q(3*n+1),q(3*n+2));

    [g_true,H_true]=cdf_grad_gx2(x,c.w,c.k,c.l,c.s,c.m,tight{:});
    tic; g_an=cdf_grad_gx2(x,c.w,c.k,c.l,c.s,c.m); t_g=toc;
    tic; [~,H_an]=cdf_grad_gx2(x,c.w,c.k,c.l,c.s,c.m); t_H=toc;

    best_g=[inf inf nan nan]; best_H=[inf inf nan nan];
    for h=steps
        tic; g=fd_grad(F,p,h); tg=toc;
        e=errs(g,g_true(:),noK(:));
        if e(1)<best_g(1), best_g=[e h tg]; end
        tic; H=fd_hess(F,p,h); tH=toc;
        e=errs(H,H_true,noK2);
        if e(1)<best_H(1), best_H=[e h tH]; end
    end

    fprintf('\n=== case %d ===\n',cs)
    report('gradient',t_g,errs(g_an(:),g_true(:),noK(:)),best_g)
    report('Hessian',t_H,errs(H_an,H_true,noK2),best_H)
end
end

function g=fd_grad(F,p,h)
P=numel(p); g=zeros(P,1);
for i=1:P
    e=zeros(1,P); e(i)=h;
    g(i)=(F(p+e)-F(p-e))/(2*h);
end
end

function H=fd_hess(F,p,h)
P=numel(p); H=zeros(P); f0=F(p);
for i=1:P
    ei=zeros(1,P); ei(i)=h;
    H(i,i)=(F(p+ei)-2*f0+F(p-ei))/h^2;
    for j=i+1:P
        ej=zeros(1,P); ej(j)=h;
        H(i,j)=(F(p+ei+ej)-F(p+ei-ej)-F(p-ei+ej)+F(p-ei-ej))/(4*h^2);
        H(j,i)=H(i,j);
    end
end
end

function e=errs(est,truth,mask)
d=abs(est-truth);
e=[max(d(:)) max(d(mask))];
end

function report(name,t_an,a,best)
fprintf('%s: analytic %.3g s, FD %.3g s (h=%.1e), speed ratio FD/an %.2g\n',...
    name,t_an,best(4),best(3),best(4)/t_an)
fprintf('  all entries : an err %.2g, FD err %.2g, acc ratio %.2g\n',a(1),best(1),best(1)/max(a(1),realmin))
fprintf('  without k   : an err %.2g, FD err %.2g, acc ratio %.2g\n',a(2),best(2),best(2)/max(a(2),realmin))
end
