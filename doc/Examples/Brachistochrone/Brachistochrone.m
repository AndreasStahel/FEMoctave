%% find a few brachistochrone curves
L = 10; u_left = 0; u_right = 1; Interval = linspace(0,L,51);

v0 = 2.0;
a = @(x,u,du)(1+du.^2).^(-1/2).*(u+v0^2).^(-1/2);
f = {@(x,u,du) 1/2*(1+du.^2).^(+1/2).*(u+v0^2).^(-3/2),
     @(x,u,du)-3/4*(1+du.^2).^(+1/2).*(u+v0^2).^(-5/2),
     @(x,u,du) 1/2*(1+du.^2).^(-1/2).*du.*(u+v0^2).^(-3/2)};

u0 = @(x)u_left+x/L*(u_right-u_left);
[x,u2,inform] = BVP1DNL(Interval,a,0,0,1,f,u_left,u_right,u0,'display','iter','maxiter',50);

v0 = 1.0;
a = @(x,u,du)(1+du.^2).^(-1/2).*(u+v0^2).^(-1/2);
f = {@(x,u,du)1/2*(1+du.^2).^(1/2).*(u+v0^2).^(-3/2),
     @(x,u,du)-3/4*(1+du.^2).^(1/2).*(u+v0^2).^(-5/2),
     @(x,u,du)1/2*(1+du.^2).^(-1/2).*du.*(u+v0^2).^(-3/2)};
[x,u1,inform] = BVP1DNL(Interval,a,0,0,1,f,u_left,u_right,u2,'display','iter','maxiter',50);

v0 = 0.75;
a = @(x,u,du)(1+du.^2).^(-1/2).*(u+v0^2).^(-1/2);
f = {@(x,u,du)1/2*(1+du.^2).^(1/2).*(u+v0^2).^(-3/2),
     @(x,u,du)-3/4*(1+du.^2).^(1/2).*(u+v0^2).^(-5/2),
     @(x,u,du)1/2*(1+du.^2).^(-1/2).*du.*(u+v0^2).^(-3/2)};
[x,u75,inform] = BVP1DNL(Interval,a,0,0,1,f,u_left,u_right,u1,'display','iter','maxiter',50);

v0 = 0.50;
a = @(x,u,du)(1+du.^2).^(-1/2).*(u+v0^2).^(-1/2);
f = {@(x,u,du)1/2*(1+du.^2).^(1/2).*(u+v0^2).^(-3/2),
     @(x,u,du)-3/4*(1+du.^2).^(1/2).*(u+v0^2).^(-5/2),
     @(x,u,du)1/2*(1+du.^2).^(-1/2).*du.*(u+v0^2).^(-3/2)};
[x,u50,inform] = BVP1DNL(Interval,a,0,0,1,f,u_left,u_right,u75,'display','iter','maxiter',50);

figure(1); plot(x,-u2,x,-u1,x,-u75,x,-u50); xlabel('x'); ylabel('-u')
           legend('v0=2.0','v0=1.0','v0=0.75','v0=0.50')
