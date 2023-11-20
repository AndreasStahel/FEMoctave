L = 60; N = 401;
Interval = linspace(0,L,N)';
w = 1; a = 1; b = 0; c = 0; d = 1; alpha = 1.0;
f = {@(x,t,u)alpha*u.*(1-u), @(x,t,u)alpha*(1-2*u)};
t0 = 0; tend = 25; steps = [30,30];
BCleft = [0,0]; BCright = 0;
u0 = @(x)0.1*exp(-3*x.^2);
[x,u_all,t] = IBVP1DNL(Interval,w,a,b,c,d,f,BCleft,BCright,u0,t0,tend,steps,'tol', 1e-8);
ind15 = find(abs(t-15)<1e-10); ind20 = find(abs(t-20)<1e-10); ind25 = find(abs(t-25)<1e-10);
figure(1); plot(x,u_all(:,ind15),x,u_all(:,ind20),x,u_all(:,ind25)); xlabel('x'); ylabel('u');
           legend('t=15','t=20','t=25','location','southwest')
figure(2); mesh(t,x,u_all); xlabel('t'); xlim([t0,tend])
           ylabel('x'); zlabel('u'); view([-70,20])
figure(3); contour(t,x,u_all); xlabel('t'); ylabel('x');
