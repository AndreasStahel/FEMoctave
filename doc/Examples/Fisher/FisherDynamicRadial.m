## -*- texinfo -*-
## @deftypefn  {} {} FisherDynamicRadial.m
##
## This is a demo file  inside the `doc/Examples/Fisher/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

L = 60; N = 401; Interval = linspace(0,L,N)';
w = @(r)r.^2; a = @(r)r.^2; b = 0; c = 0; d = @(r)r.^2;
alpha = 1.0/pi;
%%f = {@(r,t,u)alpha*u.*(1-u),@(r,t,u)alpha*(1-2*u)};
f = {@(r,t,u)alpha*sin(pi*u),@(r,t,u)alpha*pi*cos(pi*u)};
t0 = 0; tend = 25; steps = [30,60];
BCleft = [0,0]; BCright = 0;
u0 = @(r)0.1*exp(-3*r.^2);
[r,u_all,t] = IBVP1DNL(Interval,w,a,b,c,d,f,BCleft,BCright,u0,t0,tend,steps,'tol', 1e-8);
ind15 = find(abs(t-15)<1e-10); ind20 = find(abs(t-20)<1e-10); ind25 = find(abs(t-25)<1e-10);
ind5 = find(abs(t-5)<1e-10); ind10 = find(abs(t-10)<1e-10);
figure(1); plot(r,u_all(:,ind5),r,u_all(:,ind10),r,u_all(:,ind15),r,u_all(:,ind20),r,u_all(:,ind25)); xlabel('r'); ylabel('u');
           legend('t=5','t=10','t=15','t=20','t=25','location','northeast')
figure(2); mesh(t,r,u_all); xlabel('t'); xlim([t0,tend])
           ylabel('r'); zlabel('u'); view([-70,20])
figure(3); contour(t,r,u_all); xlabel('t'); ylabel('r');
