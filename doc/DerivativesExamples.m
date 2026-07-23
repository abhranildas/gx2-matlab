%% Derivatives of the generalized chi-square cdf
% Two short examples of the analytic gradient/Hessian routines added to the
% toolbox. They compute exact derivatives of the cdf (no finite differencing)
% from the characteristic-function representation.
%
% * *Example 1* differentiates in the native parameters (w, k, lambda, s, m)
%   with |cdf_grad_gx2|.
% * *Example 2* differentiates in the quadratic-boundary coefficients
%   (Q2, q1, q0) of a normal vector with |cdf_grad_norm_quad| -- the object
%   used to measure how a classification probability responds to its boundary.
%
% Make sure the toolbox is on the path:
addpath(fileparts(fileparts(mfilename('fullpath'))));

%% Example 1: gradient and Hessian in the native parameters
% Take a generalized chi-square and a point x0, and ask how its cdf value
% F(x0)=P(chi<=x0) changes as we nudge the parameters that define it.
w=[1 -5 2]; k=[1 2 3]; lambda=[2 3 7]; s=2; m=5;
x0=10;

% The gradient is a flat vector over all parameters, in the canonical order
% [w, k, lambda, s, m] (all of w, then all of k, ...); the Hessian is the
% matching square matrix.
[grad,hess]=cdf_grad_gx2(x0,w,k,lambda,s,m)

%% Taylor picture: vary one non-centrality and predict the cdf
% We isolate the derivatives with respect to lambda alone, then use the first
% component's gradient and curvature to build the second-order Taylor model of
% F(x0) as lambda(1) moves, and overlay it on the true cdf. The two touch to
% second order at delta = 0; the gap grows only as delta^3.
[g,H]=cdf_grad_gx2(x0,w,k,lambda,s,m,'wrt',{'lambda'});
gl=g(1); Hl=H(1,1);                       % d/dlambda1 and d^2/dlambda1^2
F0=gx2cdf(x0,w,k,lambda,s,m);

delta=linspace(-50,50,100);
Ftrue=arrayfun(@(d) gx2cdf(x0,w,k,lambda+[d 0 0],s,m),delta);
Ftaylor=F0+gl*delta+0.5*Hl*delta.^2;

figure; plot(lambda(1)+delta,Ftrue,'k-','LineWidth',1); hold on;
plot(lambda(1)+delta,Ftaylor,'-b','LineWidth',1);
plot(lambda(1),F0,'bo','MarkerFaceColor','b');
xlabel('\lambda_1'); ylabel('F(x_0)');
axis([-50 50 0 1])
legend('true cdf','2nd-order Taylor','location','best');
legend boxoff
title('cdf sensitivity to a non-centrality \lambda_1');

%% Example 2: sensitivity to the quadratic-boundary coefficients
% A normal vector x~N(mu,v) and a quadratic q(x)=x'*Q2*x + q1'*x + q0. The
% probability in the quadratic region up to a level x0 is F(x0)=P(q(x)<=x0);
% at x0=0 the surface q(x)=0 is a classification boundary. cdf_grad_norm_quad
% returns dF/d(Q2,q1,q0) directly, holding the normal fixed.
mu=[1;2]; v=[2 1; 1 3];
quad.q2=[1 1; 1 1]; quad.q1=[-1;0]; quad.q0=-1;
x0=0;

[grad,hess]=cdf_grad_norm_quad(x0,mu,v,quad);
disp('dF/dQ2:'); disp(grad.q2);
disp('dF/dq1:'); disp(grad.q1);
fprintf('dF/dq0: %.4f\n',grad.q0);

%% Taylor picture: vary one Q2 coefficient and predict the cdf
% We isolate the (1,1) entry of the quadratic form and use its gradient and
% curvature to build the second-order Taylor model of F(x0) as Q2(1,1) moves,
% then overlay it on the true probability. The two touch to second order at
% delta = 0; the gap grows only as delta^3.
g11=grad.q2(1,1); H11=hess.q2q2(1,1,1,1);   % dF/dQ2_11 and d^2F/dQ2_11^2
[w2,k2,l2,s2,m2]=norm_quad_to_gx2_params(mu,v,quad);
F0=gx2cdf(x0,w2,k2,l2,s2,m2);

delta=linspace(-1.5,1.5,100);
Ftrue=arrayfun(@(d) probq(mu,v,quad,d,x0),delta);
Ftaylor=F0+g11*delta+0.5*H11*delta.^2;

figure; plot(quad.q2(1,1)+delta,Ftrue,'k-','LineWidth',1); hold on;
plot(quad.q2(1,1)+delta,Ftaylor,'-b','LineWidth',1);
plot(quad.q2(1,1),F0,'bo','MarkerFaceColor','b');
xlabel('Q_2(1,1)'); ylabel('F(x_0)');
legend('true cdf','2nd-order Taylor','location','best');
legend boxoff
title('cdf sensitivity to a boundary coefficient Q_2(1,1)');

%% helper: probability with the Q2(1,1) coefficient perturbed by d
function p=probq(mu,v,quad,d,x0)
    quad.q2(1,1)=quad.q2(1,1)+d;
    [w,k,lambda,s,m]=norm_quad_to_gx2_params(mu,v,quad);
    p=gx2cdf(x0,w,k,lambda,s,m);
end
