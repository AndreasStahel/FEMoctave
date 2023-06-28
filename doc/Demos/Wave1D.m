%% test the 1D wave equation
a = 1; b = 0; c = 0; d = 1; f = 0;
w2 = 1; w1 = 0; ; BCleft = 0; BCright = [0,0];
t0 = 0; tend = 25; steps = [100,10];
interval = linspace(0,3*pi,51)';
u0 = @(x)sin(x).*(x<=pi); u1 = @(x)-cos(x).*(x<=pi);

[x,u,t] = I2BVP1D(interval,w2,w1,a,b,c,d,f,BCleft,BCright,u0,u1,t0,tend,steps);

figure(1); mesh(t,x,u); xlabel('time t'); ylabel('position x'); zlabel('u')
           xlim([min(t),max(t)]); ylim([min(x),max(x)])

figure(2); contour(t,x,u,21); xlabel('time t'); ylabel('position x');
           colorbar
