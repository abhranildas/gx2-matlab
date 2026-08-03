function bench_norm_err_grad_bd(outdir, problem_names)
% BENCH_NORM_ERR_GRAD_BD Benchmark of norm_err_grad_bd (analytic gradient/
% Hessian of total two-class classification error w.r.t. the boundary
% coefficients) on the same 23 binary-Gaussian problems as gx2-py's
% tests/bench_norm_err_bd.py, spanning D=1..5.
%
% For each problem this computes the optimal quadratic (Bayes) boundary
% (via IntClassNorm's norm_class_opt_bd), then compares three ways of
% getting the gradient/Hessian of the classification error there:
%   - a tight-tolerance analytic evaluation (norm_err_grad_bd), ground truth;
%   - the default-tolerance analytic evaluation (package defaults);
%   - finite differences (Richardson-extrapolated central differences).
% Each problem also gets the same tight/default/FD comparison for the
% gradient only (no Hessian) at a second, deliberately non-optimal
% naive-QDA boundary (the same quadratic-discriminant formula as the
% optimal boundary, but with each class's covariance replaced by its own
% diagonal), where the gradient is generically nonzero -- see
% boundary_naive_qda.
%
% No wall-clock cap: a stage runs to completion however long it needs (some
% of these problems are here specifically because one or both methods can
% be very slow -- mixed-sign, higher-dimension boundaries). Parallelism is
% across problems only (parfor); each problem's 7 stages run sequentially.
%
% Usage:
%   bench_norm_err_grad_bd()                      % all 23 problems
%   bench_norm_err_grad_bd(outdir)                 % custom output dir
%   bench_norm_err_grad_bd(outdir, {'1_D1_generic_A','4_D2_same_cov'})
%                                                   % only these problems
%                                                   % (e.g. for a quick
%                                                   % smoke test)
%
% Results: one result_<name>.json and log_<name>.log per problem, plus
% meta.json and a combined bench_norm_err_grad_bd_results.json, all written
% to outdir (default: bench_out/ next to this script).

if nargin<1 || isempty(outdir)
    outdir = fullfile(fileparts(mfilename('fullpath')), 'bench_out');
end
if ~isfolder(outdir)
    mkdir(outdir);
end

problems = get_problems();
if nargin>=2 && ~isempty(problem_names)
    keep = ismember({problems.name}, problem_names);
    problems = problems(keep);
end
n = numel(problems);

pool = gcp;
meta = struct('matlab_version', version, 'computer', computer, ...
    'n_workers', pool.NumWorkers, 'n_problems', n);
write_partial(fullfile(outdir, 'meta.json'), meta);

for i = 1:n
    p = problems(i);
    D = numel(p.mu0);
    placeholder = struct('name', p.name, 'D', D, 'P', n_params(D), ...
        'mu0', p.mu0, 'v0', p.v0, 'mu1', p.mu1, 'v1', p.v1, 'p0', p.p0, 'p1', p.p1);
    write_partial(fullfile(outdir, sprintf('result_%s.json', p.name)), placeholder);
end

parfor i = 1:n
    run_problem(problems(i), outdir);
end

combined_problems = containers.Map('KeyType', 'char', 'ValueType', 'any');
for i = 1:n
    p = problems(i);
    path = fullfile(outdir, sprintf('result_%s.json', p.name));
    if isfile(path)
        combined_problems(p.name) = jsondecode(fileread(path));
    else
        combined_problems(p.name) = struct('status', 'not started');
    end
end
combined = struct('meta', meta, 'problems', combined_problems);
write_partial(fullfile(outdir, 'bench_norm_err_grad_bd_results.json'), combined);
fprintf('Done. Combined results: %s\n', fullfile(outdir, 'bench_norm_err_grad_bd_results.json'));
end


% ---------------------------------------------------------------------------
% The problems. Transcribed verbatim from gx2-py/tests/bench_norm_err_bd.py's
% PROBLEMS list. Numbered/ordered first by dimension, then by boundary type
% (linear, elliptic, hyperbolic, parabolic) -- see gx2_derivatives.md section
% 2.2.1.
% ---------------------------------------------------------------------------
function problems = get_problems()
problems = struct('name', {}, 'mu0', {}, 'v0', {}, 'mu1', {}, 'v1', {}, 'p0', {}, 'p1', {});
k = 0;

k=k+1; problems(k) = mkprob('1_D1_generic_A', 0.0, 1.0, 3.0, 4.0);
k=k+1; problems(k) = mkprob('2_D1_generic_B', 0.0, 1.0, 1.5, 0.4);
k=k+1; problems(k) = mkprob('3_D1_near_linear', 0.0, 1.0, 2.0, 1.05);

% D=2, linear: identical covariance for both classes (Q2=0 exactly), so
% the boundary collapses to purely linear (classic LDA).
k=k+1; problems(k) = mkprob('4_D2_same_cov', [-.3;-.3], eye(2), [.3;.3], eye(2));

% D=2, elliptic: same-sign contrast -- faster same-sign code path.
k=k+1; problems(k) = mkprob('5_D2_same_sign', [-.3;-.3], eye(2), [.3;.3], diag([0.4,0.6]));

% D=2, elliptic: same covariances as problem 5, unequal priors (p0=0.2,p1=0.8).
k=k+1; problems(k) = mkprob('6_D2_unequal_prior', [-.3;-.3], eye(2), [.3;.3], diag([0.4,0.6]), 0.2, 0.8);

k=k+1; problems(k) = mkprob('7_D2_generic_axis', [-.3;-.3], eye(2), [.3;.3], diag([3.0,0.5]));
k=k+1; problems(k) = mkprob('8_D2_generic_rot', [-.3;-.3], eye(2), [.3;.3], [2.0 0.8; 0.8 1.0]);

% D=2, hyperbolic: the density-cusp case.
k=k+1; problems(k) = mkprob('9_D2_generic_crossed', [.3;-.3], diag([0.5,2.0]), [-.3;.3], diag([2.0,0.5]));

k=k+1; problems(k) = mkprob('10_D2_near_linear', [0.0;0.0], [1.0 .2; .2 1.0], [.3;.2], [1.05 .2; .2 .97]);

% D=2, hyperbolic: same mean for both classes -- q1 vanishes identically,
% boundary driven purely by the covariance difference.
k=k+1; problems(k) = mkprob('11_D2_same_mean', [0.0;0.0], eye(2), [0.0;0.0], diag([3.0,0.5]));

% D=2, parabolic: Q2 rank-deficient (one exact zero eigenvalue) but not the
% zero matrix, with q1 nonzero along that same flat axis -- a genuine
% paraboloid boundary.
k=k+1; problems(k) = mkprob('12_D2_parabolic', [-.3;0.0], eye(2), [.3;0.0], diag([1,3.0]));

% D=3 counterpart of problem 4.
k=k+1; problems(k) = mkprob('13_D3_same_cov', [-.3;-.3;-.3], eye(3), [.3;.3;.3], eye(3));

% D=3 counterpart of problem 5.
k=k+1; problems(k) = mkprob('14_D3_same_sign', [-.3;-.3;-.3], eye(3), [.3;.3;.3], diag([0.4,0.6,0.5]));

k=k+1; problems(k) = mkprob('15_D3_generic', [-.3;-.3;-.3], eye(3), [.3;.3;.3], diag([2.0,0.5,3.0]));
k=k+1; problems(k) = mkprob('16_D3_generic_crossed', [.3;0;-.3], diag([0.5,2.0,0.7]), [-.3;0;.3], diag([2.0,0.5,1.4]));
k=k+1; problems(k) = mkprob('17_D3_near_linear', [0.0;0.0;0.0], eye(3), [.3;.2;0.0], [1.05 .03 0; .03 1.02 .02; 0 .02 1.04]);

% D=3 counterpart of problem 11 (same mean, reuses problem 15's covariances).
k=k+1; problems(k) = mkprob('18_D3_same_mean', [0.0;0.0;0.0], eye(3), [0.0;0.0;0.0], diag([2.0,0.5,3.0]));

% D=3 counterpart of problem 12 (parabolic).
k=k+1; problems(k) = mkprob('19_D3_parabolic', [-.3;0.0;0.0], eye(3), [.3;0.0;0.0], diag([1,3.0,0.5]));

% D=4/D=5 same-sign vs. mixed-sign pairs, deliberately kept despite the
% mixed-sign ones being slow (see gx2_derivatives.md open item 3.2): during
% design, problem 20's covariances swapped for problem 21's mixed-sign ones
% didn't finish its analytic Hessian in 16+ minutes in the Python benchmark.
k=k+1; problems(k) = mkprob('20_D4_same_sign', -.3*ones(4,1), eye(4), .3*ones(4,1), diag([0.4,0.5,0.6,0.7]));
k=k+1; problems(k) = mkprob('21_D4_mixed_sign', -.3*ones(4,1), eye(4), .3*ones(4,1), diag([2.0,0.5,3.0,1.5]));
k=k+1; problems(k) = mkprob('22_D5_same_sign', -.3*ones(5,1), eye(5), .3*ones(5,1), diag([0.4,0.5,0.6,0.7,0.5]));
k=k+1; problems(k) = mkprob('23_D5_mixed_sign', -.3*ones(5,1), eye(5), .3*ones(5,1), diag([2.0,0.5,3.0,1.5,0.7]));
end

function p = mkprob(name, mu0, v0, mu1, v1, p0, p1)
if nargin<6, p0 = 0.5; end
if nargin<7, p1 = 0.5; end
p = struct('name', name, 'mu0', mu0(:), 'v0', v0, 'mu1', mu1(:), 'v1', v1, 'p0', p0, 'p1', p1);
end


% ---------------------------------------------------------------------------
% The naive-QDA boundary: the same quadratic-discriminant formula as
% norm_class_opt_bd, but with each class's covariance replaced by its own
% diagonal (i.e. ignoring cross-feature correlations). Unlike the Fisher/LDA
% boundary this keeps a genuinely nonzero, generically mixed-sign quadratic
% term Q2 -- so it exercises the same non-degenerate gx2 machinery
% (Ruben/Imhof, cross-component derivatives) as the optimal-boundary Hessian
% stages -- while still being deliberately non-optimal for the true (fully
% correlated) problem, so the classification error's gradient there is
% generically nonzero, unlike at the true optimum.
% ---------------------------------------------------------------------------
function quad = boundary_naive_qda(mu0, v0, mu1, v1, p0, p1)
dv0 = diag(v0); dv1 = diag(v1);
q2 = diag(0.5*(1./dv1 - 1./dv0));
q1 = mu0./dv0 - mu1./dv1;
q0 = 0.5*(sum(mu1.^2./dv1) - sum(mu0.^2./dv0)) + 0.5*(sum(log(dv1)) - sum(log(dv0))) + log(p0/p1);
quad = struct('q2', q2, 'q1', q1(:), 'q0', q0);
end


% ---------------------------------------------------------------------------
% Scalar total classification error, used as the FD objective function (and
% for the analytic stages' reported err field). Deliberately NOT
% classify_normals (IntClassNorm) -- that does a lot of unrelated work
% (boundary points, discriminability indices, sample-optimal-boundary logic)
% even with plotting off, which would make FD look artificially slow next
% to the analytic norm_err_grad_bd calls. This matches the corrected
% convention in norm_err_grad_bd.m: class 0 favored where q(x)>0.
% ---------------------------------------------------------------------------
function err = total_err(mu0, v0, mu1, v1, quad, p0, p1, varargin)
[w0,k0,l0,s0,m0] = norm_quad_to_gx2_params(mu0, v0, quad);
[w1,k1,l1,s1,m1] = norm_quad_to_gx2_params(mu1, v1, quad);
F0 = gx2cdf(0, w0, k0, l0, s0, m0, varargin{:});
F1 = gx2cdf(0, w1, k1, l1, s1, m1, varargin{:});
err = p0*F0 + p1*(1-F1);
end


% ---------------------------------------------------------------------------
% Parametrization: theta = [vech(Q2) (upper triangle incl. diagonal,
% row-major), q1, q0]. Ported directly from gx2-py/tests/bench_norm_err_bd.py
% so FD and analytic results compare component-for-component the same way.
% ---------------------------------------------------------------------------
function P = n_params(D)
P = D*(D+1)/2 + D + 1;
end

function theta = flatten_quad(quad, D)
theta = zeros(n_params(D), 1);
idx = 1;
for r = 1:D
    for c = r:D
        theta(idx) = quad.q2(r,c);
        idx = idx + 1;
    end
end
theta(idx:idx+D-1) = quad.q1(:);
theta(idx+D) = quad.q0;
end

function quad = unflatten_theta(theta, D)
q2 = zeros(D,D);
idx = 1;
for r = 1:D
    for c = r:D
        q2(r,c) = theta(idx);
        q2(c,r) = theta(idx);
        idx = idx + 1;
    end
end
q1 = theta(idx:idx+D-1);
q0 = theta(idx+D);
quad = struct('q2', q2, 'q1', q1(:), 'q0', q0);
end

function [dq2, dq1, dq0] = flatten_dir(i, D)
dq2 = zeros(D,D); dq1 = zeros(D,1); dq0 = 0.0;
idx = 1;
for r = 1:D
    for c = r:D
        if idx==i
            dq2(r,c) = 1.0; dq2(c,r) = 1.0;
        end
        idx = idx + 1;
    end
end
for j = 1:D
    if idx==i
        dq1(j) = 1.0;
    end
    idx = idx + 1;
end
if idx==i
    dq0 = 1.0;
end
end

function out = grad_flat(grad, D)
out = zeros(n_params(D), 1);
idx = 1;
for r = 1:D
    for c = r:D
        if r==c
            out(idx) = grad.q2(r,c);
        else
            out(idx) = 2*grad.q2(r,c);
        end
        idx = idx + 1;
    end
end
out(idx:idx+D-1) = grad.q1(:);
out(idx+D) = grad.q0;
end

function val = safe_contract(H, W)
% sum(H.*W) without 0*inf turning into nan -- a zero weight means "this raw
% entry isn't part of this basis direction," so it must contribute exactly
% 0 even where H itself is +-inf (as it legitimately is at a density-cusp
% boundary). W is always an exact 0/1 indicator (outer products of
% flatten_dir's 0/1 vectors), so W==0 unambiguously means "excluded."
prod = H .* W;
prod(W==0) = 0.0;
val = sum(prod(:));
end

function val = hess_dir(hess, dq2a, dq1a, dq0a, dq2b, dq1b, dq0b)
D = size(dq2a, 1);
outer_q2q2 = reshape(dq2a,[D,D,1,1]) .* reshape(dq2b,[1,1,D,D]);
val = safe_contract(hess.q2q2, outer_q2q2);
val = val + safe_contract(hess.q1q2, dq1a(:) .* reshape(dq2b,[1,D,D]));
val = val + safe_contract(hess.q1q2, dq1b(:) .* reshape(dq2a,[1,D,D]));
val = val + safe_contract(hess.q1q1, dq1a(:)*dq1b(:)');
val = val + safe_contract(hess.q0q1, dq1a(:)*dq0b) + safe_contract(hess.q0q1, dq1b(:)*dq0a);
val = val + safe_contract(hess.q0q2, dq2a*dq0b) + safe_contract(hess.q0q2, dq2b*dq0a);
val = val + safe_contract(hess.q0q0, dq0a*dq0b);
end

function H = hess_flat(hess, D)
P = n_params(D);
dq2s = cell(P,1); dq1s = cell(P,1); dq0s = cell(P,1);
for i = 1:P
    [dq2s{i}, dq1s{i}, dq0s{i}] = flatten_dir(i, D);
end
H = zeros(P,P);
for i = 1:P
    for j = i:P
        v = hess_dir(hess, dq2s{i}, dq1s{i}, dq0s{i}, dq2s{j}, dq1s{j}, dq0s{j});
        H(i,j) = v; H(j,i) = v;
    end
end
end

function v = safe_max_abs_diff(a, b)
% max(abs(a-b)) without agreeing +-inf entries turning into nan -- see
% gx2-py's _safe_max_abs_diff (same rationale, needed for problem 9's
% density-cusp Hessian).
d = abs(a - b);
agreeing_inf = isinf(a) & isinf(b) & (sign(a)==sign(b));
d(agreeing_inf) = 0.0;
v = max(d(:));
end


% ---------------------------------------------------------------------------
% Finite differences: central differences at a halving sequence of step
% sizes, refined by Richardson (Romberg-style) extrapolation -- no
% numdifftools equivalent exists in MATLAB. num_steps matches the Python
% benchmark's numdifftools settings (15 for gradient, 9 for Hessian) for
% comparable cost/precision, not a bit-identical algorithm. Progress is
% logged every 20 calls via fprintf, which the caller's diary() picks up.
% ---------------------------------------------------------------------------
function val = richardson_extrapolate(T)
n = numel(T);
R = zeros(n,n);
R(:,1) = T(:);
for j = 2:n
    for i = j:n
        R(i,j) = R(i,j-1) + (R(i,j-1) - R(i-1,j-1)) / (4^(j-1) - 1);
    end
end
val = R(n,n);
end

function [grad, n_calls] = fd_grad(fun, theta0, num_steps)
P = numel(theta0);
grad = zeros(P,1);
call_count = 0;
for i = 1:P
    h0 = 0.01*max(abs(theta0(i)), 1);
    T = zeros(num_steps,1);
    for kk = 1:num_steps
        h = h0/2^(kk-1);
        tp = theta0; tp(i) = tp(i)+h;
        tm = theta0; tm(i) = tm(i)-h;
        fp = fun(tp); call_count = call_count+1; log_fd_call(call_count, 'grad');
        fm = fun(tm); call_count = call_count+1; log_fd_call(call_count, 'grad');
        T(kk) = (fp - fm) / (2*h);
    end
    grad(i) = richardson_extrapolate(T);
end
n_calls = call_count;
end

function [hess, n_calls] = fd_hess(fun, theta0, num_steps)
P = numel(theta0);
hess = zeros(P,P);
call_count = 0;
f0 = fun(theta0); call_count = call_count+1; log_fd_call(call_count, 'hess');
for i = 1:P
    h0 = 0.01*max(abs(theta0(i)), 1);
    T = zeros(num_steps,1);
    for kk = 1:num_steps
        h = h0/2^(kk-1);
        tp = theta0; tp(i) = tp(i)+h;
        tm = theta0; tm(i) = tm(i)-h;
        fp = fun(tp); call_count = call_count+1; log_fd_call(call_count, 'hess');
        fm = fun(tm); call_count = call_count+1; log_fd_call(call_count, 'hess');
        T(kk) = (fp - 2*f0 + fm) / h^2;
    end
    hess(i,i) = richardson_extrapolate(T);
end
for i = 1:P
    for j = i+1:P
        hi0 = 0.01*max(abs(theta0(i)), 1);
        hj0 = 0.01*max(abs(theta0(j)), 1);
        T = zeros(num_steps,1);
        for kk = 1:num_steps
            hi = hi0/2^(kk-1); hj = hj0/2^(kk-1);
            tpp = theta0; tpp(i) = tpp(i)+hi; tpp(j) = tpp(j)+hj;
            tpm = theta0; tpm(i) = tpm(i)+hi; tpm(j) = tpm(j)-hj;
            tmp = theta0; tmp(i) = tmp(i)-hi; tmp(j) = tmp(j)+hj;
            tmm = theta0; tmm(i) = tmm(i)-hi; tmm(j) = tmm(j)-hj;
            fpp = fun(tpp); call_count = call_count+1; log_fd_call(call_count, 'hess');
            fpm = fun(tpm); call_count = call_count+1; log_fd_call(call_count, 'hess');
            fmp = fun(tmp); call_count = call_count+1; log_fd_call(call_count, 'hess');
            fmm = fun(tmm); call_count = call_count+1; log_fd_call(call_count, 'hess');
            T(kk) = (fpp - fpm - fmp + fmm) / (4*hi*hj);
        end
        v = richardson_extrapolate(T);
        hess(i,j) = v; hess(j,i) = v;
    end
end
n_calls = call_count;
end

function log_fd_call(call_count, kind)
if mod(call_count, 20)==0
    fprintf('[%s]   FD %s call %d\n', datestr(now,'HH:MM:SS'), kind, call_count);
end
end


% ---------------------------------------------------------------------------
% Analytic stages: total_err for the reported err field, norm_err_grad_bd
% (direct call, natural mu0/mu1 order -- see the sign-convention fix in
% norm_err_grad_bd.m) for gradient/Hessian.
% ---------------------------------------------------------------------------
function nv = tol_to_namevalue(tol)
nv = {};
if isfield(tol,'AbsTol'), nv = [nv, {'AbsTol', tol.AbsTol}]; end
if isfield(tol,'RelTol'), nv = [nv, {'RelTol', tol.RelTol}]; end
if isfield(tol,'n_ruben'), nv = [nv, {'n_ruben', tol.n_ruben}]; end
end

function [err_v, grad_flat_v, hess_flat_v] = analytic_stage(mu0, v0, mu1, v1, quad, p0, p1, tol, D)
nv = tol_to_namevalue(tol);
err_v = total_err(mu0, v0, mu1, v1, quad, p0, p1, nv{:});
[grad_raw, hess_raw] = norm_err_grad_bd(mu0, v0, mu1, v1, quad, 'p0', p0, 'p1', p1, nv{:});
grad_flat_v = grad_flat(grad_raw, D);
hess_flat_v = hess_flat(hess_raw, D);
end

function [err_v, grad_flat_v] = analytic_stage_grad_only(mu0, v0, mu1, v1, quad, p0, p1, tol, D)
nv = tol_to_namevalue(tol);
err_v = total_err(mu0, v0, mu1, v1, quad, p0, p1, nv{:});
grad_raw = norm_err_grad_bd(mu0, v0, mu1, v1, quad, 'p0', p0, 'p1', p1, nv{:});
grad_flat_v = grad_flat(grad_raw, D);
end


% ---------------------------------------------------------------------------
% Per-problem worker: boundary setup, then all 7 stages sequentially (4
% optimal-boundary, 3 naive-QDA-boundary). Writes result_<name>.json to disk
% after every stage (atomic tmp+rename) so a killed run leaves usable
% partial results. diary() captures every fprintf (and any MATLAB warning)
% printed during this problem's run into log_<name>.log -- one mechanism
% for both "verbose log files" and "duplicate command-line output."
% ---------------------------------------------------------------------------
function run_problem(problem, outdir)
% here = <repo>/gx2-matlab/tests: gx2-matlab itself (for norm_quad_to_gx2_params,
% gx2cdf, norm_err_grad_bd, ...) is one level up; IntClassNorm (a sibling of
% the gx2 repo itself, for norm_class_opt_bd) is three levels up.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..', '..', 'IntClassNorm'));

name = problem.name;
log_path = fullfile(outdir, sprintf('log_%s.log', name));
result_path = fullfile(outdir, sprintf('result_%s.json', name));

diary(log_path);
diary on;

mu0 = problem.mu0; v0 = problem.v0; mu1 = problem.mu1; v1 = problem.v1;
p0 = problem.p0; p1 = problem.p1;
D = numel(mu0);
P = n_params(D);

result = struct('name', name, 'D', D, 'P', P, 'mu0', mu0, 'v0', v0, ...
    'mu1', mu1, 'v1', v1, 'p0', p0, 'p1', p1);
write_partial(result_path, result);
fprintf('[%s] starting, D=%d, P=%d\n', datestr(now,'HH:MM:SS'), D, P);

try
    quad = norm_class_opt_bd([mu0,v0], [mu1,v1], 'prior_1', p0);
    theta0 = flatten_quad(quad, D);
    quad_naive_qda = boundary_naive_qda(mu0, v0, mu1, v1, p0, p1);
    theta0_naive_qda = flatten_quad(quad_naive_qda, D);
    result.quad = quad;
    result.quad_naive_qda = quad_naive_qda;
    write_partial(result_path, result);
    fprintf('[%s] boundaries computed (optimal + naive-QDA)\n', datestr(now,'HH:MM:SS'));
catch ME
    result.error = getReport(ME);
    write_partial(result_path, result);
    fprintf('[%s] PROBLEM_FAILED (boundary setup):\n%s\n', datestr(now,'HH:MM:SS'), result.error);
    diary off;
    return;
end

tight = struct('AbsTol', 1e-12, 'RelTol', 1e-10, 'n_ruben', 5000);
stages = {'ground_truth','default_tol','fd_grad','fd_hess', ...
    'naive_qda_ground_truth','naive_qda_default_tol','naive_qda_fd_grad'};

for si = 1:numel(stages)
    stage = stages{si};
    tic_id = tic;
    status_str = 'ok';
    value = [];
    try
        switch stage
            case 'ground_truth'
                [ev,gv,hv] = analytic_stage(mu0,v0,mu1,v1,quad,p0,p1,tight,D);
                value = struct('err',ev,'grad',gv,'hess',hv);
            case 'default_tol'
                [ev,gv,hv] = analytic_stage(mu0,v0,mu1,v1,quad,p0,p1,struct(),D);
                value = struct('err',ev,'grad',gv,'hess',hv);
            case 'fd_grad'
                [g,nc] = fd_grad(@(th) total_err(mu0,v0,mu1,v1,unflatten_theta(th,D),p0,p1), theta0, 15);
                value = struct('grad',g,'n_calls',nc);
            case 'fd_hess'
                [h,nc] = fd_hess(@(th) total_err(mu0,v0,mu1,v1,unflatten_theta(th,D),p0,p1), theta0, 9);
                value = struct('hess',h,'n_calls',nc);
            case 'naive_qda_ground_truth'
                [ev,gv] = analytic_stage_grad_only(mu0,v0,mu1,v1,quad_naive_qda,p0,p1,tight,D);
                value = struct('err',ev,'grad',gv);
            case 'naive_qda_default_tol'
                [ev,gv] = analytic_stage_grad_only(mu0,v0,mu1,v1,quad_naive_qda,p0,p1,struct(),D);
                value = struct('err',ev,'grad',gv);
            case 'naive_qda_fd_grad'
                [g,nc] = fd_grad(@(th) total_err(mu0,v0,mu1,v1,unflatten_theta(th,D),p0,p1), theta0_naive_qda, 15);
                value = struct('grad',g,'n_calls',nc);
        end
    catch ME
        status_str = 'error';
        value = getReport(ME);
    end
    elapsed = toc(tic_id);
    result = store_stage_result(result, stage, status_str, value, elapsed, D);
    write_partial(result_path, result);
    fprintf('[%s] %s: %s in %.2fs\n', datestr(now,'HH:MM:SS'), stage, status_str, elapsed);
end

result = finalize_problem(result);
write_partial(result_path, result);
fprintf('[%s] PROBLEM_DONE\n', datestr(now,'HH:MM:SS'));
diary off;
end


function result = store_stage_result(result, stage, status, value, elapsed, D) %#ok<INUSD>
switch stage
    case 'ground_truth'
        if strcmp(status,'ok')
            result.ground_truth = struct('status',status,'time_s',elapsed,'err',value.err, ...
                'grad',value.grad,'hess',value.hess,'tol','AbsTol=1e-12, RelTol=1e-10, n_ruben=5000');
        else
            result.ground_truth = struct('status',status,'time_s',elapsed,'detail',value);
        end
    case 'default_tol'
        if strcmp(status,'ok')
            result.default_tol = struct('status',status,'time_s',elapsed,'err',value.err, ...
                'grad',value.grad,'hess',value.hess,'tol','AbsTol=1e-10, RelTol=1e-6 (package defaults)');
        else
            result.default_tol = struct('status',status,'time_s',elapsed,'detail',value);
        end
    case 'fd_grad'
        if strcmp(status,'ok')
            result.fd_grad = struct('status',status,'time_s',elapsed,'grad',value.grad, ...
                'n_calls',value.n_calls,'fd_opts','Richardson-extrapolated central differences, num_steps=15');
        else
            result.fd_grad = struct('status',status,'time_s',elapsed,'detail',value);
        end
    case 'fd_hess'
        if strcmp(status,'ok')
            result.fd_hess = struct('status',status,'time_s',elapsed,'hess',value.hess, ...
                'n_calls',value.n_calls,'fd_opts','Richardson-extrapolated central differences, num_steps=9');
        else
            result.fd_hess = struct('status',status,'time_s',elapsed,'detail',value);
        end
    case {'naive_qda_ground_truth','naive_qda_default_tol'}
        if strcmp(stage,'naive_qda_ground_truth')
            fkey = 'ground_truth'; tol = 'AbsTol=1e-12, RelTol=1e-10, n_ruben=5000';
        else
            fkey = 'default_tol'; tol = 'AbsTol=1e-10, RelTol=1e-6 (package defaults)';
        end
        if ~isfield(result,'naive_qda'), result.naive_qda = struct(); end
        if strcmp(status,'ok')
            result.naive_qda.(fkey) = struct('status',status,'time_s',elapsed,'err',value.err,'grad',value.grad,'tol',tol);
        else
            result.naive_qda.(fkey) = struct('status',status,'time_s',elapsed,'detail',value);
        end
    case 'naive_qda_fd_grad'
        if ~isfield(result,'naive_qda'), result.naive_qda = struct(); end
        if strcmp(status,'ok')
            result.naive_qda.fd_grad = struct('status',status,'time_s',elapsed,'grad',value.grad, ...
                'n_calls',value.n_calls,'fd_opts','Richardson-extrapolated central differences, num_steps=15');
        else
            result.naive_qda.fd_grad = struct('status',status,'time_s',elapsed,'detail',value);
        end
end
end


function result = finalize_problem(result)
grad_gt = []; hess_gt = [];
if isfield(result,'ground_truth') && strcmp(result.ground_truth.status,'ok')
    grad_gt = result.ground_truth.grad;
    hess_gt = result.ground_truth.hess;
end

analytic_err_grad = []; analytic_err_hess = [];
if isfield(result,'default_tol') && strcmp(result.default_tol.status,'ok')
    if ~isempty(grad_gt)
        analytic_err_grad = safe_max_abs_diff(result.default_tol.grad, grad_gt);
        result.default_tol.analytic_err_grad = analytic_err_grad;
    end
    if ~isempty(hess_gt)
        analytic_err_hess = safe_max_abs_diff(result.default_tol.hess, hess_gt);
        result.default_tol.analytic_err_hess = analytic_err_hess;
    end
end

fd_err_grad = [];
if isfield(result,'fd_grad') && strcmp(result.fd_grad.status,'ok') && ~isempty(grad_gt)
    fd_err_grad = max(abs(result.fd_grad.grad(:) - grad_gt(:)));
    result.fd_grad.err = fd_err_grad;
end

fd_err_hess = [];
if isfield(result,'fd_hess') && strcmp(result.fd_hess.status,'ok') && ~isempty(hess_gt)
    fd_err_hess = max(abs(result.fd_hess.hess(:) - hess_gt(:)));
    result.fd_hess.err = fd_err_hess;
end

t_default = []; if isfield(result,'default_tol') && strcmp(result.default_tol.status,'ok'), t_default = result.default_tol.time_s; end
t_fd_grad = []; if isfield(result,'fd_grad') && strcmp(result.fd_grad.status,'ok'), t_fd_grad = result.fd_grad.time_s; end
t_fd_hess = []; if isfield(result,'fd_hess') && strcmp(result.fd_hess.status,'ok'), t_fd_hess = result.fd_hess.time_s; end

total_time = 0;
for f = {'ground_truth','default_tol','fd_grad','fd_hess'}
    if isfield(result,f{1}) && isfield(result.(f{1}),'time_s')
        total_time = total_time + result.(f{1}).time_s;
    end
end
if isfield(result,'naive_qda')
    fn = fieldnames(result.naive_qda);
    for fi = 1:numel(fn)
        if isfield(result.naive_qda.(fn{fi}),'time_s')
            total_time = total_time + result.naive_qda.(fn{fi}).time_s;
        end
    end
end

summary = struct('total_stage_time_s', total_time);
summary.rel_speed_grad = ratio_or_empty(t_fd_grad, t_default);
summary.rel_speed_hess = ratio_or_empty(t_fd_hess, t_default);
summary.rel_acc_grad = ratio_or_empty(fd_err_grad, analytic_err_grad);
summary.rel_acc_hess = ratio_or_empty(fd_err_hess, analytic_err_hess);
result.summary = summary;

if isfield(result,'naive_qda')
    naive_qda = result.naive_qda;
    qgrad_gt = [];
    if isfield(naive_qda,'ground_truth') && strcmp(naive_qda.ground_truth.status,'ok')
        qgrad_gt = naive_qda.ground_truth.grad;
    end
    naive_qda_analytic_err_grad = [];
    if isfield(naive_qda,'default_tol') && strcmp(naive_qda.default_tol.status,'ok') && ~isempty(qgrad_gt)
        naive_qda_analytic_err_grad = max(abs(naive_qda.default_tol.grad(:) - qgrad_gt(:)));
        result.naive_qda.default_tol.analytic_err_grad = naive_qda_analytic_err_grad;
    end
    naive_qda_fd_err_grad = [];
    if isfield(naive_qda,'fd_grad') && strcmp(naive_qda.fd_grad.status,'ok') && ~isempty(qgrad_gt)
        naive_qda_fd_err_grad = max(abs(naive_qda.fd_grad.grad(:) - qgrad_gt(:)));
        result.naive_qda.fd_grad.err = naive_qda_fd_err_grad;
    end
    t_naive_qda_default = []; if isfield(naive_qda,'default_tol') && strcmp(naive_qda.default_tol.status,'ok'), t_naive_qda_default = naive_qda.default_tol.time_s; end
    t_naive_qda_fd = []; if isfield(naive_qda,'fd_grad') && strcmp(naive_qda.fd_grad.status,'ok'), t_naive_qda_fd = naive_qda.fd_grad.time_s; end
    qsummary = struct();
    qsummary.rel_speed_grad = ratio_or_empty(t_naive_qda_fd, t_naive_qda_default);
    qsummary.rel_acc_grad = ratio_or_empty(naive_qda_fd_err_grad, naive_qda_analytic_err_grad);
    result.naive_qda.summary = qsummary;
end
end

function r = ratio_or_empty(num, den)
if isempty(num) || isempty(den) || den==0
    r = [];
else
    r = num/den;
end
end


function write_partial(path, result)
tmp = [path '.tmp'];
txt = jsonencode(result, 'PrettyPrint', true);
fid = fopen(tmp, 'w');
fwrite(fid, txt);
fclose(fid);
movefile(tmp, path, 'f');
end
