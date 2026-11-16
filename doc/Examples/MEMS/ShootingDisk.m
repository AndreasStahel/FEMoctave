R = 0.5; u0 = 0.01;

function res = IVP(lambda,u0)
  R = 0.5;
  options = odeset('InitialStep', 1e-4,'MaxStep',1e-3);
  ODE = @(r,u) [u(2);-u(2)/r-lambda./(1-u(1)).^2];
  [r,u] = ode45(ODE,[0.0001,R],[u0;0]);%,options);
  res = u(end,1);
endfunction


ToSolve = @(lambda)IVP(lambda,u0);
lambda = fsolve(ToSolve,0.01)
N = 50;
u0_List = linspace(0.01,0.9,N); lambda_List = zeros(size(u0_List));

for jj = 1:length(u0_List)
  ToSolve = @(lambda)IVP(lambda,u0_List(jj));
  [lambda,u_value] = fsolve(ToSolve,u0_List(jj));
  lambda_List(jj) = lambda;
endfor
figure(5); plot(lambda_List,u0_List);
           xlabel('\lambda'); ylabel('u(0)');

% plot an upper solution
lambda = lambda_List(N); u0 = u0_List(N);
ODE = @(r,u) [u(2);-u(2)/r-lambda./(1-u(1)).^2];
[r,u] = ode45(ODE,[0.0001,R],[u0;0]);
figure(6); plot(r,u(:,1)); xlabel('r'); ylabel('u(x)')
u = u(:,1);
save UpperSolutionDisk.mat r u lambda u0
% plot a lower solution
N2 = round(N/2);
lambda = lambda_List(N2); u0 = u0_List(N2);
ODE = @(r,u) [u(2);-u(2)/r-lambda./(1-u(1)).^2];
[r,u] = ode45(ODE,[0.0001,R],[u0;0]);
figure(7); plot(r,u(:,1)); xlabel('r'); ylabel('u(x)')
