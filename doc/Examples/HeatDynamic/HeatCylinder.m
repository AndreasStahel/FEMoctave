R = 3; BCleft = [0,0]; BCright = [0,0];
f = @(r,t)sin(1*t)*(r<R/2);
u0 = 0;  t0 = 0; t_end = 6*pi; steps = [121,10];  interval = linspace(0,R,21);
r = @(r)r; solver = 'CN';
[r,u,t] = IBVP1D(interval,r,r,0,0,r,f,BCleft,BCright,u0,t0,t_end,steps,'solver',solver);

figure(1); mesh(t,r,u)
           xlabel('time t'); ylabel('radius r'); zlabel('temperature u')
figure(2); contour(t,r,u,[-0.5:0.1:+0.5])
           xlabel('time t'); ylabel('radius r');

figure(3); plot(t,u(1,:),t,u(end,:))
           xlabel('time t'); ylabel('temperature u'); xlim([0,max(t)])
           legend('at center r=0','at outer edge r=R','location','south')

