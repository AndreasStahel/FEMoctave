## -*- texinfo -*-
## @deftypefn  {} {} SoundWaveSpherical1D.m
##
## This is a demo file  inside the `doc/Examples/SoundWave/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

N = 2*60; R = 3; interval = linspace(0,R,N)';
f_r2 = @(r)r.^2;
u0 = @(r)(1+cos(10*r)).*(r<pi/10);   u1 = 0;
w2 = f_r2; w1 = 0; a = f_r2; b = 0; c = 0; d = 1; f = 0;
BCleft = [0,0]; BCright = [0,0];
t0 = 0; tend = 2.5; steps = [100,10];
[r,u,t] = I2BVP1D(interval,w2,w1,a,b,c,d,f,BCleft,BCright,u0,u1,t0,tend,steps);

figure(1); mesh(t,r,u); xlabel('time t'); ylabel('radius r'); zlabel('u')
           xlim([min(t),max(t)]); ylim([min(r),max(r)])

[t05,t05_ind] = find(abs(t-0.5)<1e-5); [t1,t1_ind] = find(abs(t-1)<1e-5);
[t15,t15_ind] = find(abs(t-1.5)<1e-5); [t2,t2_ind] = find(abs(t-2)<1e-5);
[t25,t25_ind] = find(abs(t-2.5)<1e-5);

figure(2); plot(r,u(:,1),r,u(:,t05_ind),r,u(:,t1_ind),r,u(:,t15_ind),r,u(:,t2_ind),r,u(:,t25_ind));
           ylabel('amplitude u'); ylim([-0.5,0.755]); xlabel('radius r');
           legend('t=0.0','t=0.5','t=1.0','t=1.5','t=2.0','t=2.5','location','northeast')
