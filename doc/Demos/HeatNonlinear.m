pkg load femoctave
L = 1; N = 41;  Interval = linspace(0,L,N)';
w = 1; a = 0.1 ; b = 0; c = 0; d = 1;
f = {@(x,t,u)x.^3+sin(u),@(x,t,u)cos(u)};
t0 = 0; tend = 10; steps = [30,50];
BCleft = 0; BCright = 0;
u0 = 0;
[x,u_all,t] = IBVP1DNL(Interval,w,a,b,c,d,f,BCleft,BCright,u0,t0,tend,steps);
figure(1); mesh(t,x,u_all); xlabel('t'); xlim([t0,tend])
           ylabel('x'); zlabel('u'); view([30,30])
figure(2); contour(t,x,u_all); xlabel('t'); ylabel('x');

[x,u] = BVP1DNL(Interval,a,b,c,d,{@(x,u)x.^3+sin(u),@(x,u)cos(u)},BCleft,BCright,u0);
figure(3); plot(x,u,x,u_all(:,end)); xlabel('x'); ylabel('u');
           legend('static','dynamic', 'location','south')
