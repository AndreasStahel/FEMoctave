## -*- texinfo -*-
## @deftypefn  {} {} WaveExplicitTest.m
##
## This is a demo file inside the `doc/Demos/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

%% test the stability 1D wave equation
if 1
 solver = 'explicit';
else
 solver = 'implicit';
endif
solver = 'expl'

if 0
  steps = [20,5];  %% unstable for explicit solver
else
  steps = [21,5];  %% stable for both solvers
endif

a = 1; b = 0; c = 0; d = 1; f = 0; w2 = 1; w1 = 0; BCleft = 0; BCright = [0,0];
t0 = 0; tend = 5; interval = linspace(0,3*pi,51)';
u0 = @(x)sin(x).*(x<=pi); u1 = @(x)-cos(x).*(x<=pi);
[x,u,t] = I2BVP1D(interval,w2,w1,a,b,c,d,f,BCleft,BCright,u0,u1,t0,tend,steps,'solver',solver);

figure(11); clf; mesh(t,x,u); xlabel('time t'); ylabel('position x'); zlabel('u')
            xlim([min(t),max(t)]); ylim([min(x),max(x)])
figure(12); clf; contour(t,x,u,21); xlabel('time t'); ylabel('position x');
figure(13); plot(x,u(:,end)); xlabel('x'); ylabel('u(x) at t_{end}'); xlim([min(x),max(x)])
