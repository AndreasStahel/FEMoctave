N = 51; interval = linspace(-1,1,N)';
c = 1.176501939901833; %% or use the solver
opt.TolFun = 1e-15; opt.TolX = 1e-15;  c = fsolve(@(c)1+cos(c)-c^2,1,opt);
CASE = 3;  %% select between successive substitution (1) and Newton's method (2)
           %% or the function BVP1DNL()  (3)
switch CASE
case 1  %% successive substitution
  [x,u] = BVP1D(interval,1,0,0,1,-1,0,0);
  u_exact = log(c^2./(1+cos(c*x)));
  figure(1); plot(x,u); xlabel('x'); ylabel('u_1')
  for jj = 1:4;
    [x,u] = BVP1D(interval,1,0,0,1,-exp(u),0,0);
    figure(2); plot(x,u,x,u_exact); xlabel('x'); ylabel('u')
               legend('FEM','exact', 'location','north')
    figure(3); plot(x,u-u_exact); xlabel('x'); ylabel('u')
               legend('FEM-exact')
    pause(0.2)
  endfor
case 2  %% Newton's method
  [x,u] = BVP1D(interval,1,0,0,1,-1,0,0);
  u_exact = log(c^2./(1+cos(c*x)));
  figure(1); plot(x,u); xlabel('x'); ylabel('u_1')
  xGauss = FEM1DGaussPoints(x);
  for jj = 1:4
    [du,ddu] = FEM1DEvaluateDu(x,u);
    RHS = + ddu - exp(u);
    uGauss = pwquadinterp(x,u,xGauss); %% evaluate u at Gauss points
    u_coeff = exp(uGauss);
    [x,phi] = BVP1D(interval,1,0,u_coeff,1,RHS,0,0);
    disp(sprintf('max(abs(phi)) = %g, max(abs(RHS)) = %g',max(abs(phi)), max(abs(RHS))))
    u = u + phi;
    figure(2); plot(x,u,x,u_exact); xlabel('x'); ylabel('u')
               legend('FEM','exact','location','north')
    figure(3); plot(x,u-u_exact); xlabel('x'); ylabel('u')
               legend('FEM-exact')
    pause(0.2)
  endfor
case 3 %% use BVP1DNL
  f = {@(x,u)-exp(u) , @(x,u)-exp(u)};
  [x,u,inform] = BVP1DNL(interval,1,0,0,1,f,0,0,0,'Tol',1e-5,'display','iter');
  u_exact = log(c^2./(1+cos(c*x)));
  figure(2); plot(x,u,x,u_exact); xlabel('x'); ylabel('u')
             legend('FEM','exact','location','north')
  figure(3); plot(x,u-u_exact); xlabel('x'); ylabel('u')
             legend('FEM-exact')
  inform
endswitch
RMS_difference = norm(u-u_exact)/sqrt(length(u))

