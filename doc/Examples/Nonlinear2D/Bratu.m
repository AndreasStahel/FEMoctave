## -*- texinfo -*-
## @deftypefn  {} {} Bratu.m
##
## This is a example file inside the `doc/Examples/Nonlinear2D/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

C = 1.0; C_crit = 3.513830719;
N = 201;
Interval = linspace(0,1,N)';
f = {@(x,u)C*exp(u),@(x,u)C*exp(u)};
u0 = @(x)0.8*sin(pi*x); %% BVP1DNL, one solution, lower branch
[x,u_low] = BVP1DNL(Interval,1,0,0,1,f,0,0,u0,'display','iter');
u_max = [max(u_low),pwquadinterp(x,u_low,0.5)]

u0 = @(x)4*sin(pi*x); %% BVP1DNL, one solution, upper branch
[x,u_upp] = BVP1DNL(Interval,1,0,0,1,f,0,0,u0,'display','iter');
u_max = [max(u_upp),pwquadinterp(x,u_upp,0.5)]

u0 = @(x)1.2*sin(pi*x);
f = {@(x,u)C_crit*exp(u),@(x,u)C_crit*exp(u)};
[x,u_crit] = BVP1DNL(Interval,1,0,0,1,f,0,0,u0,'display','iter','maxiter',20);
u_max = [max(u_crit),pwquadinterp(x,u_crit,0.5)]
figure(1); plot(x,[u_low,u_upp,u_crit]); xlabel('x'); ylabel('u')
           legend('C=1.0 lower','C=1.0 upper','C=C_{crit}')
