## -*- texinfo -*-
## @deftypefn  {} {} Reflection.m
##
## This is a demo file  inside the `doc/Examples/WaveAnimation/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

interval = linspace(0,14,301)';
a = @(x) 1+1*(x>4);
b = 0; c = 0; d = 1; f = 0;
w2 = 1; w1 = 0; ; BCleft = [0,0]; BCright = [0,0];
t0 = 0; tend = 9; steps = [90,100];
u0 = @(x)sin(pi*x).*(x<=1); u1 = @(x)-pi*cos(pi*x).*(x<=1);

[x,u,t] = I2BVP1D(interval,w2,w1,a,b,c,d,f,BCleft,BCright,u0,u1,t0,tend,steps);

figure(1); mesh(t,x,u); xlabel('time t'); ylabel('position x'); zlabel('u')
           xlim([min(t),max(t)]); ylim([min(x),max(x)]); view([-20,10])

t_ind = [1 21 41 61 81];
figure(2); H = waterfall(x,t(t_ind),u(:,t_ind)'); view([-177,35])
           xlabel('position x'); ylabel('t'); zlabel('u')
           set(H, 'edgecolor', [0,0.0,1])
           set(H, 'linewidth',3)

figure(3); contour(t,x,u,[-0.1,0.1]); xlabel('time t'); ylabel('position x');

