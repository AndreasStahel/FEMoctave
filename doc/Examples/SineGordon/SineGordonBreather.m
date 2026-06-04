## -*- texinfo -*-
## @deftypefn  {} {} SineGordonBreather.m
##
## This is a demo file inside the `doc/Examples/SineGordon/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

omega = 0.5; L = 10; Interval = linspace(0,L,501); Tend = 4.5*pi/omega;

f  = @(x,t,u)-sin(u);
v0 = @(x)4*sqrt(1-omega^2)./cosh(sqrt(1-omega^2)*x);
[x,u,t] = I2BVP1DNL(Interval,1,0,1,0,0,1,f,[0,0],0,0,v0,0,Tend,[61,100] ,'Solver','implicit');

figure(1); mesh(t,x,u); view(-60,15)
           xlabel('time t'); ylabel('position x')
figure(2); contour(t,x,u,21)
           xlabel('time t'); ylabel('position x')

u_max = @(x)4*atan(sqrt(1-omega^2)./(omega*cosh(sqrt(1-omega^2)*x)));
figure(3); plot(x,u(:,end),x,u_max(x))
           xlabel('position x'); ylabel('u'); legend('FEM','exact')
figure(4); plot(x,u_max(x)-u(:,end))
           xlabel('position x'); ylabel("u_{exact}-u_{FEM}");
