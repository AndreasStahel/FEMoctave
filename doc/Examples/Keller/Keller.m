N = 51;
flag = 1;  %% select between successive substition (1) and Newton's method (2)
switch flag
case 1  %% successive substitution
  x = linspace(-1,1,N);
  [xn,u] = BVP1D(x,1,0,0,1,-1,0,0);
  c = 1.1765019; u_exact = log(c^2./(1+cos(c*xn)));
  figure(1); plot(xn,u); xlabel('x'); ylabel('u_1')
  for jj = 1:4;
    [xn,u] = BVP1D(x,1,0,0,1,-exp(u),0,0);
    figure(2); plot(xn,u,xn,u_exact); xlabel('x'); ylabel('u')
               legend('FEM','exact', 'location','north')
    figure(3); plot(xn,u-u_exact); xlabel('x'); ylabel('u')
               legend('FEM-exact')
    pause(1)
  endfor
case 2  %% Newton's method
  x = linspace(-1,1,N);
  [xn,u] = BVP1D(x,1,0,0,1,-1,0,0);
  c = 1.1765019; u_exact = log(c^2./(1+cos(c*xn)));
  figure(1); plot(xn,u); xlabel('x'); ylabel('u_1')
  xGauss = FEM1DGaussPoints(xn);
  for jj = 1:4
    [du,ddu] = FEM1DEvaluateDu(xn,u);
    RHS = ddu - exp(u);
    uGauss = pwquadinterp(xn,u,xGauss); %% evaluate u at Gauss points
    u_coeff = exp(uGauss);
    [xn,phi] = BVP1D(x,1,0,u_coeff,1,RHS,0,0);
    maxPHI_RHS= [max(abs(phi)) max(abs(RHS))]
    u = u + phi;
    figure(2); plot(xn,u,xn,u_exact); xlabel('x'); ylabel('u')
               legend('FEM','exact', 'location','north')
    figure(3); plot(xn,u-u_exact); xlabel('x'); ylabel('u')
               legend('FEM-exact')
    pause(1)
  endfor
endswitch

