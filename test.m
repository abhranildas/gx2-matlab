
%% debug IntClassNorm
mu_1=[4; 5]; v_1=[2 1; 1 1];
mu_2=mu_1; v_2=3*[2 -1; -1 1];

results_gx2=classify_normals([mu_1,v_1],[mu_2,v_2],'method','gx2')

%% pdf of vector function
mu=[-1; -1]; v=[1 0.5; 0.5 2];
fun=@(x,y) x.*sin(y) - y.*cos(x);

fun_grad=@(x) [sin(x(2))+x(2)*sin(x(1)); x(1)*cos(x(2))-cos(x(1))];

x=linspace(-20,20,50);
% p=norm_fun_cdf(x,mu,v,fun);
% plot(x,p)
tic
f_diff=norm_fun_pdf(x,mu,v,fun,'pdf_method','diff','dx',0.1,'AbsTol',0,'RelTol',1e-5);
time_diff=toc

tic
f_ray=norm_fun_pdf(x,mu,v,fun,'pdf_method','ray','fun_grad',fun_grad,'AbsTol',0,'RelTol',1e-5);
time_ray=toc
figure; hold on
plot(x,f_diff)
plot(x,f_ray)

%% new gx2pdf_ray
w=[-2 1 1];
k=[2 1 2];
l=[7 5 3];
m=5;
s=3;

x=linspace(-60,40,1e2);
tic
f=gx2pdf_ray(x,w,k,l,s,m,'mc_samples',1e5);
toc

plot(x,f)

% TODO near the peak of the pdf it's noisy with MC,
% and with quadrature it uses a lot of rays. Why?

%% new gx2pdf_ray, grid integral, 2d

w=[-2 1];
k=[1 1];
l=[7 5];
m=5;
s=0;

x=linspace(-60,40,1e2);
tic
f=gx2pdf_ray(x,w,k,l,s,m);
toc
% f2=gx2pdf_imhof(x,w,k,l,s,m); %TODO Imhof bad in this and next

figure; hold on
plot(x,f)
plot(x,f2)

%% test gx2pdf with 'method'

w=[-2 1 3];
k=[1 1 7];
l=[7 5 0];
m=0;
s=10;

x=linspace(-60,40,1e2);

% tic
% f_auto=gx2pdf(x,w,k,l,s,m);
% toc

tic
f_imhof=gx2pdf(x,w,k,l,s,m,'method','imhof');
toc

tic
f_ray=gx2pdf(x,w,k,l,s,m,'method','ray','n_rays',1e5,'vpa',false);
toc

% TODO why is ray output wiggly around -20


figure; hold on
% plot(x,f_auto)
plot(x,f_imhof,'-k')
plot(x,f_ray,'-b')
% plot(x,f_conv)

%% gx2inv with log prob

% non-elliptic gx2

w=[1 -10 2];
k=[1 2 3];
l=[2 3 7];
s=5;
m=10;

x=gx2inv(-1e3,w,k,l,s,m,'method','ray')
p=gx2cdf(x,w,k,l,s,m,'method','ray')

x=gx2inv(-1e3,w,k,l,s,m,'upper','method','ray')
p=gx2cdf(x,w,k,l,s,m,'upper','method','ray')

% elliptic gx2
w=[4 5 1];
k=[1 2 3];
l=[5 6 0];
s=0;
m=-50;

x=gx2inv(-1e3,w,k,l,s,m,'method','ellipse','x_scale','log')
p=gx2cdf(x,w,k,l,s,m,'method','ellipse','x_scale','log')

x=gx2inv(-1e3,w,k,l,s,m,'upper','method','tail')
p=gx2cdf(x,w,k,l,s,m,'upper','method','ray')


%% test tail pdf and cdf approximation

w=[-2 -3 4 1];
k=[4 2 6 1];
l=[10 20 5 0];
s=20;
m=0;

x=linspace(0,2e3,100);

% upper tail

% pdf
log10_f_ray=gx2pdf(x,w,k,l,s,m,'method','ray','n_rays',1e3);
log10_f_ray(log10_f_ray>0)=log10(log10_f_ray(log10_f_ray>0));
f_imhof=gx2pdf(x,w,k,l,s,m,'method','imhof','abstol',1e-100,'reltol',1e-100);
f_pearson=gx2pdf(x,w,k,l,s,m,'method','pearson');
log10_f_tail=gx2pdf(x,w,k,l,s,m,'upper','method','tail');
log10_f_tail(log10_f_tail>0)=log10(log10_f_tail(log10_f_tail>0));

figure; hold on
plot(x,log10(f_imhof),'-ok')
plot(x,log10_f_ray,'-ob')
plot(x,log10(f_pearson),'-og')
plot(x,log10_f_tail,'-or')
yline(log10(realmin))

% cdf
log10_p_ray=gx2cdf(x,w,k,l,s,m,'upper','method','ray','n_rays',1e7);
log10_p_ray(log10_p_ray>0)=log10(log10_p_ray(log10_p_ray>0));
p_imhof=gx2cdf(x,w,k,l,s,m,'upper','method','imhof','abstol',1e-100,'reltol',1e-100);
p_pearson=gx2cdf(x,w,k,l,s,m,'upper','method','pearson');
log10_p_tail=gx2cdf(x,w,k,l,s,m,'upper','method','tail');
log10_p_tail(log10_p_tail>0)=log10(log10_p_tail(log10_p_tail>0));

figure; hold on
plot(x,log10(p_imhof),'ok')
plot(x,log10_p_ray,'-b')
plot(x,log10(p_pearson),'-g')
plot(x,log10_p_tail,'or')
yline(log10(realmin))

% lower tail

x=-linspace(0,8e3,20);

% pdf
log10_f_ray=gx2pdf(x,w,k,l,s,m,'method','ray','n_rays',1e7);
log10_f_ray(log10_f_ray>0)=log10(log10_f_ray(log10_f_ray>0));
f_imhof=gx2pdf(x,w,k,l,s,m,'method','imhof','abstol',1e-200,'reltol',1e-200);
f_pearson=gx2pdf(x,w,k,l,s,m,'method','pearson');
log10_f_tail=gx2pdf(x,w,k,l,s,m,'lower','method','tail');
log10_f_tail(log10_f_tail>0)=log10(log10_f_tail(log10_f_tail>0));

figure; hold on
plot(x,log10(f_imhof),'-ok')
plot(x,log10_f_ray,'-ob')
plot(x,log10(f_pearson),'-g')
plot(x,log10_f_tail,'-or')
yline(log10(realmin))

% cdf
log10_p_ray=gx2cdf(x,w,k,l,s,m,'method','ray','n_rays',1e7);
log10_p_ray(log10_p_ray>0)=log10(log10_p_ray(log10_p_ray>0));
p_imhof=gx2cdf(x,w,k,l,s,m,'method','imhof','abstol',1e-200,'reltol',1e-200);
p_pearson=gx2cdf(x,w,k,l,s,m,'method','pearson');
log10_p_tail=gx2cdf(x,w,k,l,s,m,'lower','method','tail');
log10_p_tail(log10_p_tail>0)=log10(log10_p_tail(log10_p_tail>0));

figure; hold on
plot(x,log10(p_imhof),'ok')
plot(x,log10_p_ray,'-b')
plot(x,log10(p_pearson),'-g')
plot(x,log10_p_tail,'-r')
yline(log10(realmin))

%% log sum exp approximation: infinite tail

w=[-2 -3 4 1];
k=[4 2 6 1];
l=[10 20 5 0];
s=20;
m=0;

x=linspace(0,1e4,20);

% cdf
p_ray_log=gx2cdf(x,w,k,l,s,m,'upper','method','ray','n_rays',1e3,'gpu_batch',0);


figure; hold on
plot(x(p_ray_log>0),log10(p_ray_log(p_ray_log>0)),'-ob','MarkerFaceColor','w')
plot(x(p_ray_log<0),p_ray_log(p_ray_log<0),'-ob','MarkerFaceColor','b')

% pdf
f_ray_log=gx2pdf(x,w,k,l,s,m,'method','ray','n_rays',1e3,'gpu_batch',0);


figure; hold on
plot(x(p_ray_log>0),log10(p_ray_log(p_ray_log>0)),'-ob','MarkerFaceColor','w')
plot(x(p_ray_log<0),p_ray_log(p_ray_log<0),'-ob','MarkerFaceColor','b')

%% log sum exp approximation: finite tail

w=[2 3 4 1];
k=[4 2 6 1];
l=[0 0 0 0];
s=0;
m=0;

x=logspace(-70,-40,20);

% cdf
tic
[p_ray_log,p_err]=gx2cdf(x,w,k,l,s,m,'method','ray','n_rays',1e2);
toc
% tic
% p_ray_vpa=gx2cdf(x,w,k,l,s,m,'method','ray','precision','vpa','n_rays',10);
% toc

figure; hold on
plot(x(p_ray_log>0),log10(p_ray_log(p_ray_log>0)),'-ob','MarkerFaceColor','w')
plot(x(p_ray_log<0),p_ray_log(p_ray_log<0),'-ob','MarkerFaceColor','b')
set(gca,'xscale','log')
% plot(x,double(log10(cell2sym(p_ray_vpa))),'-or')