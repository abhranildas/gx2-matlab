<p align="center">
  <img src="https://raw.githubusercontent.com/abhranildas/gx2-matlab/main/gx2_icon.png" alt="gx2" width="260">
</p>

# Generalized chi-square distribution [![View Generalized chi-square distribution on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/85028-generalized-chi-square-distribution)
Matlab toolbox that computes the statistics, characteristic function, pdf, cdf, inverse cdf,
random numbers, and exact gradients/Hessians of the cdf, of the **generalized chi-square distribution**. There is also a [python package](https://pypi.org/project/gx2/) version of this.

A generalized chi-square variable is a weighted sum of independent non-central
chi-square variables plus a normal variable — equivalently, the quadratic form
of a normal random vector. It is parametrized by:

| parameter | meaning |
|-----------|---------|
| `w`       | weights of the non-central chi-square terms |
| `k`       | their degrees of freedom |
| `lambda` | their non-centralities |
| `s`       | scale (standard deviation) of the added normal term |
| `m`       | constant offset |

## Author and citation

Abhranil Das, Center for Perceptual Systems, The University of Texas at Austin.
Bugs / comments / questions / suggestions to abhranil.das@utexas.edu.

If you use this code, please cite:
 - [A method to integrate and classify normal distributions](https://doi.org/10.1167/jov.21.10.1)
 - [New methods to compute the generalized chi-square distribution](https://www.tandfonline.com/doi/abs/10.1080/00949655.2025.2501401)

## Installation
Within Matlab's Home tab, select Add-Ons > Get Add-Ons > Search for 'Generalized chi-square distribution' and install.

## Public functions

| function | purpose |
|----------|---------|
| `gx2stat(w, k, lambda, s, m)` | mean and variance |
| `gx2char(t, w, k, lambda, s, m)` | characteristic function |
| `gx2rnd(w, k, lambda, s, m, sz, method=)` | random numbers |
| `gx2cdf(x, w, k, lambda, s, m, side=, method=, ...)` | cdf |
| `gx2pdf(x, w, k, lambda, s, m, side=, method=, ...)` | pdf |
| `gx2inv(p, w, k, lambda, s, m, side=, method=, ...)` | inverse cdf |
| `gx2_to_norm_quad_params(w, k, lambda, s, m)` | gx2 → quadratic-form coefficients of a standard normal |
| `norm_quad_to_gx2_params(mu, v, quad, merge=)` | quadratic form of a normal → gx2 parameters |
| `cdf_grad_gx2(x, w, k, lambda, s, m, wrt=, ...)` | exact gradient (and, as a 2nd output, Hessian) of the cdf wrt the native parameters `w, k, lambda, s, m` |
| `cdf_grad_norm_quad(x, mu, v, quad, wrt=, ...)` | exact gradient (and, as a 2nd output, Hessian) of the cdf wrt the quadratic boundary coefficients `q2, q1, q0` |

For full documentation of any function, type in Matlab, e.g.:

    doc gx2_to_norm_quad_params
    doc norm_quad_to_gx2_params
    doc gx2stat
    doc gx2rnd
    doc gx2char
    doc gx2cdf
    doc gx2pdf
    doc gx2inv
    doc cdf_grad_gx2
    doc cdf_grad_norm_quad

## Computation methods for `cdf` / `pdf`

`method='auto'` (default) picks a good method for the given parameters. You can
also force one:

| method | notes |
|--------|-------|
| `'imhof'`   | Imhof–Davies numerical integration (`precision='basic'` or `'vpa'`) |
| `'ray'`     | ray-trace method (`precision='basic'`, `'log'` or `'vpa'`; tune with `n_rays`, `force_mc`) |
| `'ifft'`    | inverse-FFT method; `x='full'` returns the cdf/pdf over a spanning grid |
| `'ruben'`   | Ruben's series — requires all `w` the same sign and `s=0` |
| `'tail'`    | infinite-tail approximation |
| `'pearson'` | Pearson's 3-moment approximation |
| `'ellipse'` | ellipse approximation near a finite tail — requires all `w` the same sign and `s=0` |

## Examples
After installation, begin with the Getting Started live script with interactive examples, or, at any time, go to Matlab Home tab > Add-Ons > Manage Add-Ons > click the three dots next to this toolbox > View Getting Started Guide.

The following are the worked examples from it.

<!-- BEGIN GENERATED: getting-started (do not edit by hand; regenerate with `python scripts/build_getting_started.py`) -->

### Calculate mean and variance

```matlab
% gx2 parameters
w=[1 -10 2];
k=[1 2 3];
lambda=[2 3 7];
s=5;
m=10;

[mu,v]=gx2stat(w,k,lambda,s,m)
```
```
mu =
-17

v =
1771
```

### Generate random samples

```matlab
r=gx2rnd(w,k,lambda,s,m,[1 1e5]);
figure;
histogram(r,'EdgeColor','none')
```
![Plot output 1](https://raw.githubusercontent.com/abhranildas/gx2-matlab/main/getting-started/generate-random-samples-1.png)

### Compute PDF, CDF and inverse CDF with default methods

```matlab
x=[10 25];
f=gx2pdf(x,w,k,lambda,s,m)
```
```
f =
    0.0121    0.0088
```
```matlab
p=gx2cdf(x,w,k,lambda,s,m)
```
```
p =
    0.7150    0.8790
```
```matlab
% find the median by using the inverse CDF function:
x_med=gx2inv(.5,w,k,lambda,s,m)
```
```
x_med =
-8.7657
```
```matlab
% Compute quantiles for cdf values of 1e-3 and 1e-2, by supplying their log10 values:
x_q=gx2inv([-3 -2],w,k,lambda,s,m)
```
```
Warning: Imhof method output(s) too close to limit to compute exactly, so clipping. Check the flag output, and try stricter tolerances.

x_q =
 -219.2553 -149.2618
```
```matlab
% verify that cdf values here are indeed 1e-3 and 1e-2
p=gx2cdf(x_q,w,k,lambda,s,m)
```
```
p =
    0.0010    0.0100
```
```matlab
% Compute quantiles for complementary cdf values of 1e-3 and 1e-2, by supplying their log10 values:
x_q=gx2inv([-3 -2],w,k,lambda,s,m,'upper')
```
```
x_q =
   69.4899   51.0338
```
```matlab
% verify that ccdf values here are indeed 1e-3 and 1e-2
p=gx2cdf(x_q,w,k,lambda,s,m,'upper')
```
```
p =
    0.0010    0.0100
```
```matlab
% compute the PDF over most of the span of the distribution.
% with the 'full' argument, the span x is computed automatically.
[f,~,x]=gx2pdf('full',w,k,lambda,s,m);

% now compare the sampled histogram with the computed PDF
figure; hold on
plot(x,f)
histogram(r,'normalization','pdf','displaystyle','stairs')
xline(x_med,'-',{'median'},'labelorientation','aligned') % mark the computed median
xlim([-250 100])
```
![compute the PDF over most of the span of the distribution.](https://raw.githubusercontent.com/abhranildas/gx2-matlab/main/getting-started/compute-pdf-cdf-and-inverse-cdf-with-default-methods-2.png)
```matlab
% compute CDF over most of the span of the distribution.
% the 'full' argument uses the IFFT method, good for quick rough plots,
% but less accurate (esp. for CDF) than some other methods
[p,~,x]=gx2cdf('full',w,k,lambda,s,m);

% now compare the sampled histogram with the computed CDF
figure; hold on
plot(x,p)
histogram(r,'normalization','cdf','displaystyle','stairs')
% mark the computed median, and verify that it sits at 0.5 on the vertical axis:
xline(x_med,'-',{'median'},'labelorientation','aligned')
yline(0.5)
xlim([-200 100])
```
![compute CDF over most of the span of the distribution.](https://raw.githubusercontent.com/abhranildas/gx2-matlab/main/getting-started/compute-pdf-cdf-and-inverse-cdf-with-default-methods-3.png)

### Compute CDF, PDF and inverse CDF with each exact method and its settings

#### A non-elliptic distribution

```matlab
w=[-2 -5 2];
k=[2 1 3];
lambda=[0 4 4];
s=3;
m=-20;

% first find the quantile points at 0.1% in each tail
x_bounds=gx2inv([0.001 0.999],w,k,lambda,s,m)
```
```
Warning: Imhof method output(s) too close to limit to compute exactly, so clipping. Check the flag output, and try stricter tolerances.

x_bounds =
 -142.8108   24.8079
```
```matlab
% now compute within this range
x=linspace(x_bounds(1),x_bounds(2),50);

% compute CDF
p_ifft=gx2cdf(x,w,k,lambda,s,m,'method','ifft');
p_imhof=gx2cdf(x,w,k,lambda,s,m,'method','imhof');
p_ray=gx2cdf(x,w,k,lambda,s,m,'method','ray','n_rays',1e4);
figure; hold on
plot(x,p_ifft,'-k')
plot(x,p_imhof,'.b')
plot(x,p_ray,'or')
legend('IFFT','Imhof','ray')
```
![now compute within this range](https://raw.githubusercontent.com/abhranildas/gx2-matlab/main/getting-started/a-non-elliptic-distribution-4.png)
```matlab
% compute PDF
f_ifft=gx2pdf(x,w,k,lambda,s,m,'method','ifft');
f_imhof=gx2pdf(x,w,k,lambda,s,m,'method','imhof');
f_ray=gx2pdf(x,w,k,lambda,s,m,'method','ray','n_rays',1e6);

figure; hold on
plot(x,f_ifft,'-k')
plot(x,f_imhof,'.b')
plot(x,f_ray,'or')
legend('IFFT','Imhof','ray')
```
![compute PDF](https://raw.githubusercontent.com/abhranildas/gx2-matlab/main/getting-started/a-non-elliptic-distribution-5.png)
```matlab
% Compute quantiles for tiny cdf values of 1e-1000 and 1e-2000, by supplying
% their log10 values. Use a forward cdf method that can get down to such tiny values.
% Here we use the infinite-tail approximation.
x_q=gx2inv([-1e3 -2e3],w,k,lambda,s,m,'method','tail')
```
```
Warning: Some output values are too small for double precision. Returning their log10 values, which are negative.

x_q =
10^4 ×
   -2.4365   -4.7950
```
```matlab
% now verify using an exact cdf method that cdf values here are indeed 1e-1000 and 1e-2000:
p=gx2cdf(x_q,w,k,lambda,s,m,'method','ray','n_rays',1e7)
```
```
Warning: Some output values are too small for double precision. Returning their log10 values, which are negative.

p =
10^3 ×
   -1.0013   -2.0043
```
```matlab
% now do the same for the upper tail:
x_q=gx2inv([-1e3 -2e3],w,k,lambda,s,m,'upper','method','tail')
```
```
Warning: Some output values are too small for double precision. Returning their log10 values, which are negative.

x_q =
10^4 ×
    0.9724    1.9159
```
```matlab
p=gx2cdf(x_q,w,k,lambda,s,m,'upper','method','ray','n_rays',1e7)
```
```
Warning: Some output values are too small for double precision. Returning their log10 values, which are negative.

p =
10^3 ×
   -1.0006   -2.0012
```

#### An elliptic distribution

Here we can use Ruben's method too.

```matlab
w=[3 4 5];
k=[1 2 3];
lambda=[2 3 7];
s=0;
m=-100;

% first find the quantile points at 0.1% in each tail
x_bounds=gx2inv([0.001 0.999],w,k,lambda,s,m)
```
```
x_bounds =
  -90.5258  122.3200
```
```matlab
% now compute within this range
x=linspace(x_bounds(1),x_bounds(2),50);

% compute CDF
p_ifft=gx2cdf(x,w,k,lambda,s,m,'method','ifft');
p_imhof=gx2cdf(x,w,k,lambda,s,m,'method','imhof');
p_ray=gx2cdf(x,w,k,lambda,s,m,'method','ray','n_rays',1e4);
p_ruben=gx2cdf(x,w,k,lambda,s,m,'method','ruben');

figure; hold on
plot(x,p_ifft,'-k')
plot(x,p_imhof,'.b')
plot(x,p_ray,'or')
plot(x,p_ruben,'og','MarkerSize',8)
legend('IFFT','Imhof','ray', 'Ruben')
```
![now compute within this range](https://raw.githubusercontent.com/abhranildas/gx2-matlab/main/getting-started/an-elliptic-distribution-6.png)
```matlab
% compute PDF
f_ifft=gx2pdf(x,w,k,lambda,s,m,'method','ifft');
f_imhof=gx2pdf(x,w,k,lambda,s,m,'method','imhof');
f_ray=gx2pdf(x,w,k,lambda,s,m,'method','ray','n_rays',1e6);
f_ruben=gx2pdf(x,w,k,lambda,s,m,'method','ruben');

figure; hold on
plot(x,f_ifft,'-k')
plot(x,f_imhof,'.b')
plot(x,f_ray,'or')
plot(x,f_ruben,'og','MarkerSize',8)
legend('IFFT','Imhof','ray', 'Ruben')
```
![compute PDF](https://raw.githubusercontent.com/abhranildas/gx2-matlab/main/getting-started/an-elliptic-distribution-7.png)
```matlab
% Compute quantiles for tiny cdf values of 1e-1000 and 1e-2000, by supplying
% their log10 values. Use a forward cdf method that can get down to such tiny values.
% Here we use the ellipse approximation, with 'x_scale', 'log', which allows to specify
% log10 values of x measured from the finite tail m.
x_q=gx2inv([-1e3 -2e3],w,k,lambda,s,m,'method','ellipse','x_scale','log')
```
```
x_q =
 -331.2746 -664.6080
```
```matlab
% this means that the computed quantiles are 1e-331 and 1e-664 above m

% now verify using the forward cdf method that cdf values here are indeed 1e-1000 and 1e-2000:
p=gx2cdf(x_q,w,k,lambda,s,m,'method','ellipse','x_scale','log')
```
```
p =
10^3 ×
   -1.0000   -2.0000
```

### Compute CDF and PDF in the far tails, using some tail approximation methods too

Ray, tail and Imhof methods are best for infinite tails.

#### Compute CDF in an infinite lower tail

```matlab
w=[1 2 -3 -4];
k=[6 5 4 3];
lambda=[5 10 0 0];
s=10;
m=-50;

x=linspace(-500,200,40);

p_ifft=gx2cdf(x,w,k,lambda,s,m,'method','ifft','span',1e7,'n_grid',1e7);
p_imhof=gx2cdf(x,w,k,lambda,s,m,'method','imhof','abstol',0,'reltol',1e-10);
```
```
Warning: Imhof method output(s) too close to limit to compute exactly, so clipping. Check the flag output, and try stricter tolerances.
```
```matlab
p_ray=gx2cdf(x,w,k,lambda,s,m,'method','ray','n_rays',1e6);
p_pearson=gx2cdf(x,w,k,lambda,s,m,'method','pearson'); % pearson sucks

% tail approximation for lower tail. Mentioning 'lower' is needed here.
% For output values that are too small for double precision, it returns
% their log10 values, which are negative.
p_tail=gx2cdf(x,w,k,lambda,s,m,'lower','method','tail');
% convert all output values to their log10
p_tail(p_tail>0)=log10(p_tail(p_tail>0));

figure; hold on
plot(x,log10(p_ifft),'-k')
plot(x,log10(p_ray),'or')
plot(x,p_tail,'-g')
plot(x,log10(p_pearson),'.c','MarkerSize',15)
plot(x,log10(p_imhof),'.b')

axis([-5e2 200 -30 0])

legend('IFFT','ray','tail', 'pearson','Imhof')
ylabel('$\log_{10} p$','Interpreter','latex')
```
![tail approximation for lower tail. Mentioning 'lower' is needed here.](https://raw.githubusercontent.com/abhranildas/gx2-matlab/main/getting-started/compute-cdf-in-an-infinite-lower-tail-8.png)

#### Compute PDF in an infinite upper tail

```matlab
x=linspace(0,500,40);

f_ifft=gx2pdf(x,w,k,lambda,s,m,'method','ifft','span',1e7,'n_grid',1e7);
f_imhof=gx2pdf(x,w,k,lambda,s,m,'method','imhof','abstol',0,'reltol',1e-1);
```
```
Warning: Imhof method output(s) too close to limit to compute exactly, so clipping. Check the flag output, and try stricter tolerances.
```
```matlab
f_ray=gx2pdf(x,w,k,lambda,s,m,'method','ray','n_rays',1e6);
f_pearson=gx2pdf(x,w,k,lambda,s,m,'method','pearson');

% tail approximation for upper tail. Mentioning 'upper' is needed here.
f_tail=gx2pdf(x,w,k,lambda,s,m,'upper','method','tail');

figure; hold on
plot(x,log10(f_ifft),'-k')
plot(x,log10(f_ray),'or')
plot(x,log10(f_tail),'-g')
plot(x,log10(f_pearson),'.c','MarkerSize',15)
plot(x,log10(f_imhof),'.b')

axis([0 500 -30 0])

legend('IFFT','ray','tail', 'pearson','Imhof')
ylabel('$\log_{10} f$','Interpreter','latex')
```
![tail approximation for upper tail. Mentioning 'upper' is needed here.](https://raw.githubusercontent.com/abhranildas/gx2-matlab/main/getting-started/compute-pdf-in-an-infinite-upper-tail-9.png)

#### Compute CDF in a finite lower tail

Ruben and ellipse methods are best for finite tails.

```matlab
w=[1 2 3 4];
k=[6 5 4 3];
lambda=[5 10 0 0];
s=0;
m=0;

x=logspace(-2,2,40);

p_ifft=gx2cdf(x,w,k,lambda,s,m,'method','ifft','span',1e7,'n_grid',1e7);
p_imhof=gx2cdf(x,w,k,lambda,s,m,'method','imhof','abstol',0,'reltol',1e-10);
```
```
Warning: Imhof method output(s) too close to limit to compute exactly, so clipping. Check the flag output, and try stricter tolerances.
```
```matlab
p_ruben=gx2cdf(x,w,k,lambda,s,m,'method','ruben');
p_ray=gx2cdf(x,w,k,lambda,s,m,'method','ray','n_rays',1e5);
p_pearson=gx2cdf(x,w,k,lambda,s,m,'method','pearson');
p_ellipse=gx2cdf(x,w,k,lambda,s,m,'method','ellipse');

figure; hold on
plot(x,log10(p_ifft),'-k')
plot(x,log10(p_ray),'or','MarkerSize',8)
plot(x,log10(p_ellipse),'-g')
plot(x,log10(p_pearson),'.c','MarkerSize',15)
plot(x,log10(p_imhof),'.b')
plot(x,log10(p_ruben),'om','MarkerSize',4)

set(gca,'xscale','log')
legend('IFFT','ray','ellipse', 'pearson','Imhof','Ruben','Location', 'southeast')
ylabel('$\log_{10} p$','Interpreter','latex')
```
![Plot output 10](https://raw.githubusercontent.com/abhranildas/gx2-matlab/main/getting-started/compute-cdf-in-a-finite-lower-tail-10.png)

### Compute characteristic function

```matlab
t=linspace(-1,1,1e3);
phi=gx2char(t,w,k,lambda,s,m);
figure; plot(phi,'-o')
```
![Plot output 11](https://raw.githubusercontent.com/abhranildas/gx2-matlab/main/getting-started/compute-characteristic-function-11.png)

### Distribution of quadratic form of a normal variable

Normal parameters:

```matlab
mu=[5;6]; % mean
v=[2 1; 1 3]; % covariance matrix
```

Sample normal random vectors:

```matlab
x=mvnrnd(mu,v,1e5)';
figure; plot(x(1,:),x(2,:),'.')
```
![Plot output 12](https://raw.githubusercontent.com/abhranildas/gx2-matlab/main/getting-started/distribution-of-quadratic-form-of-a-normal-variable-12.png)

Quadratic form $q(\mathbf{x})=(x_1+x_2)^2-x_1-1$ = [x1;x2]'*[1 1; 1 1]*[x1;x2] + [-1;0]'*[x1;x2] -1

```matlab
quad.q2=[1 1; 1 1];
quad.q1=[-1;0];
quad.q0=-1;
```

Compute the quadratic form q for the sample of normal vectors:

```matlab
q=dot(x,quad.q2*x)+quad.q1'*x+quad.q0;
```

Get generalized chi-square parameters corresponding to this quadratic form:

```matlab
[w,k,lambda,s,m]=norm_quad_to_gx2_params(mu,v,quad)
```
```
w =
7.0000

k =
1

lambda =
16.6188

s =
0.8452

m =
-1.3316
```

Compare the sampled and calculated distributions of q:

```matlab
[f,~,x]=gx2pdf('full',w,k,lambda,s,m);
plot(x,f); hold on
histogram(q,'normalization','pdf','displaystyle','stairs')
xlim([0 400])
```
![Plot output 13](https://raw.githubusercontent.com/abhranildas/gx2-matlab/main/getting-started/distribution-of-quadratic-form-of-a-normal-variable-13.png)

Compare the sampled and calculated means and variances:

```matlab
[mu_q,v_q]=gx2stat(w,k,lambda,s,m);
[mu_q mean(q)]
```
```
ans =
  122.0000  121.8776
```
```matlab
[v_q var(q)]
```
```
ans =
10^3 ×
    3.3560    3.3457
```

Compare the sampled and calculated probabilities $p(q(\mathbf{x})<50)$:

```matlab
mean(q<50)
```
```
ans =
0.0858
```
```matlab
gx2cdf(50,w,k,lambda,s,m)
```
```
ans =
0.0856
```

Find a canonical quadratic form of a standard multinormal corresponding to these generalized chi-square parameters:

```matlab
quad=gx2_to_norm_quad_params(w,k,lambda,s,m)
```
```
quad =
struct with fields:
    q2: [2×2 double]
    q1: [2×1 double]
    q0: 115
```

<!-- END GENERATED: getting-started -->
