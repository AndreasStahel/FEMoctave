%% test partial successive substitution with shortest connection
interval = linspace(0,1,21)';
a = @(x,u,du) 1./sqrt(1+du.^2);
CASE = 1;
switch CASE
  case 1   %%
    u0 = @(x)2-cos(2*pi*x);
    [x,u,inform] = BVP1DNL(interval,a,0,0,1,0,1,2,u0,
                          'Tol',1e-4,'Display','iter');
    figure(1); plot(x,u); xlabel('x'); ylabel('u')
    inform
  case 2
    alpha = 0.5;
    u0 = @(x)1+x-0.3*(1-cos(2*pi*x));
    [x,u,inform] = BVP1DNL(interval,a,0,0,1,0,1,[alpha,0],u0,
                          'tol',1e-6,'MaxIter',50,'tol',1e-4,'Display','iter');
    figure(1); plot(x,u0(x),x,u); xlabel('x'); ylabel('u')
               legend('u_0','u','location','northwest')
    inform
    du = FEM1DEvaluateDu(x,u);
    u_end_du_end = [u(end) du(end) du(end)*a(0,0,du(end))]
endswitch

