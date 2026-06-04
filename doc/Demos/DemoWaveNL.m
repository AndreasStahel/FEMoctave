## -*- texinfo -*-
## @deftypefn  {} {} DemoWaveNL.m
##
## This is a demo file  inside the `doc/Demos/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

L = 10; Interval = linspace(0,L,101); Tend = 50; D = 1;
f  = @(x,t,u)0.5*(exp(-x/L)-0.5*sin(u));
v0 = @(x)0.1*x.*sin(pi*x/L);

[x,u,t] = I2BVP1DNL(Interval,1,D,1,0,0,1,f,0,0,0,v0,0,Tend,[51,100] ,'Solver','implicit');
figure(1); surf(t,x,u); view(-60,15)
           xlabel('time t'); ylabel('position x')
figure(2); contour(t,x,u,21)
           xlabel('time t'); ylabel('position x')

u0 = @(x)sin(pi*x/L);
f_stab  = {@(x,u)0.5*(exp(-x/L)-0.5*sin(u)), @(x,u)-0.5*0.5*cos(u)};
[x,u_stat] = BVP1DNL(Interval,1,0,0,1,f_stab,0,0,u0);
figure(3); plot(x,u(:,end),x,u_stat)
           xlabel('position x'); ylabel('u')
figure(4); plot(x,u(:,end)-u_stat)
           xlabel('position x'); ylabel('u-u_{stat}')

