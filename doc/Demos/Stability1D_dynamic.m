BCleft = [0,0]; BCright = 1; xMax = 1;
interval = linspace(0,xMax,15);
t0 = 0 ; tend = 2; steps = [10,1];
w = 1; a = 1; b = 0; c = 0; d = 1; f = 0;
solver = 'implicit';
%% solver = 'explicit';
%% solver = 'CN';
solver = 'RK';
[x,u,t] = IBVP1D(interval,w,a,b,c,d,f,BCleft,BCright,u0,t0,tend,steps,'solver',solver);

figure(1); clf; mesh(t,x,u); ylim([0,xMax])
           xlabel('t'); ylabel('x'); zlabel('u')
u_at_0 = u(1,end)

u_ones = ones(size(x));
figure(2); plot(x,u_ones,x,u(:,end)); xlim([0,xMax])
           xlabel('x'); ylabel('u')
           legend('1', 'FEM','location','west')

