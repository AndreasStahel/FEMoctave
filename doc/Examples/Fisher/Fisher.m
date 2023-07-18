c = 2.0; rho = 1;
M = 4; N = 21; interval = linspace(0,M,N)';
BCleft = 0.5; BCright = 1;
[x,u0] = BVP1D(interval,1,0,0,1,0,BCleft,BCright);
f = {@(x,u)rho*(-u.*(1-u)), @(x,u)rho*(-1+2*u)};
[x,u,inform] = BVP1DNL(interval,1,-c,0,1,f,BCleft,BCright,u0,
                       'Display','iter','tol',1e-8);
inform
figure(1); plot(x,u0,x,u); xlabel('x'); ylabel('u');
           legend('u_0','u','location','northwest')
figure(2); plot(x,u,'b',-x, 1-u,'b'); xlabel('x'); ylabel('u');

