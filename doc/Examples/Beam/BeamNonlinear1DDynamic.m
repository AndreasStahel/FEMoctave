%% solve the nonlinear beam problem as dynamic problem
L = 3; EI = 1; F2 = 2.0;
Interval = linspace(0,L,100)';
T0 = 0; Tend = 4; steps = [10,1];
f = {@(x,t,u)F2*cos(u), @(x,t,u)-F2*sin(u)};
u0 = 0;                             %% the naive initial guess
%%u0 = @(s,t)pi/2*(1-(s/L-1).^2);     %% with a better initial guess
%%u0 = @(s,t)-5/4*pi*(1-(s/L-1).^2);  %% to aim for a different solution
[s,u_all,t] = IBVP1DNL(Interval,1,EI,0,0,1,f,0,[0,0],u0,T0,Tend,steps);
figure(1); mesh(t,s,u_all); view([-50,20])
           xlabel('time t'); ylabel('arclength s'); zlabel('angle \alpha(s)')

u = u_all(:,end);
x = cumtrapz(s,cos(u)); y = cumtrapz(s,sin(u));
figure(2); plot(x,y); xlabel('x'); ylabel('y');
