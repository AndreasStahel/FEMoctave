R = 3; BCleft = [0,0]; BCright = 0;  f = 1;
u0 = 0;  t0 = 0; t_end = 3;
steps = [10,10];  interval = linspace(0,R,11);
r_square = @(r) r.^2;
solver = 'CN';
[r,u,t] = IBVP1D(interval,r_square,r_square,0,0,r_square,f,BCleft,BCright,...
                 u0,t0,t_end,steps, 'solver',solver);

figure(1); plot(r,u(:,end))
           xlabel('radius r'); ylabel('temperature u at t=t_{end}')
figure(2); plot(t,u(1,:))
           xlabel('time t'); ylabel('temperature u at r=R')
figure(3); mesh(t,r,u)
           xlabel('time t'); ylabel('radius r'); zlabel('temperature u')
figure(4); contour(t,r,u,[0.25:0.25:1.5])
           xlabel('time t'); ylabel('radius r');

