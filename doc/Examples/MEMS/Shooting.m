
u0 = 0.01;

function res = IVP(lambda,u0)
  ODE = @(x,u) [u(2);-lambda./(1-u(1)).^2];
  [x,u] = ode45(ODE,[0.5,1],[u0;0]);
  res = u(end,1);
endfunction

ToSolve = @(lambda)IVP(lambda,u0);
lambda = fsolve(ToSolve,0.01);
N = 100;
u0_List = linspace(0.01,0.9,N); lambda_List = zeros(size(u0_List));
for jj = 1:length(u0_List)
  ToSolve = @(lambda)IVP(lambda,u0_List(jj));
  lambda = fsolve(ToSolve,u0_List(jj));
  lambda_List(jj) = lambda;
endfor
figure(5); plot(lambda_List,u0_List);
           xlabel('\lambda'); ylabel('u(0.5)'); xlim([0,1.5])

% plot an upper solution
lambda = lambda_List(N); u0 = u0_List(N);
ODE = @(x,u) [u(2);-lambda./(1-u(1)).^2];
[x,u] = ode45(ODE,[0.5,1],[u0;0]);
figure(7); plot(x,u(:,1)); xlabel('x'); ylabel('u(x)')
x_all = linspace(0,1,1001); u_all = interp1(x,u(:,1),0.5+abs(x_all-0.5));
figure(6); plot(x_all,u_all); xlabel('x'); ylabel('u(x)')
save UpperSolution.mat x_all u_all lambda u0

printing = 0
if printing
  figure(5); print -dpdfcrop MEMSparameters.pdf
  figure(6); print -dpdfcrop MEMSUpperSolution.pdf
endif

