pkg load femoctave

CASE = 2;

switch CASE
  case 1  %% construct the lower branch
N = 51; interval = linspace(0,1,N)';
lambda_List = [[0.2:0.05:1.4]];
lambda = lambda_List(1);
f = {@(x,u)lambda./(1-u).^2 , @(x,u)+lambda*2./(1-u).^3};
[x,u,inform] = BVP1DNL(interval,1,0,0,1,f,0,0,0,'Tol',1e-5,'display','iter');
figure(1); clf; plot(x,u); xlabel('x'); ylabel('u')
hold on
u_mesh = zeros(2*N-1,length(lambda_List));
u_mesh(:,1) = u;
for jj = 2:length(lambda_List)
  lambda = lambda_List(jj)
  f = {@(x,u)lambda./(1-u).^2, @(x,u)+lambda*2./(1-u).^3};
  [x,u,inform] = BVP1DNL(interval,1,0,0,1,f,0,0,u,'Tol',1e-5,'display','iter');
  figure(1); plot(x,u); xlabel('x'); ylabel('u')
  u_mesh(:,jj) = u;
endfor
hold off
figure(2); waterfall(x,lambda_List,u_mesh'); title('lower branch');
           xlabel('x'); ylabel('\lambda'); zlabel('u(x)'); view(40,45)
case 2  %% construct the upper branch
load UpperSolution.mat;
N = 51; interval = linspace(0,1,N)';
x = sort([interval; interval(1:end-1) + diff(interval)/2]);
u = interp1(x_all,u_all,x); u0 = u; lambda0 = lambda;
f = {@(x,u)lambda./(1-u).^2, @(x,u)+lambda*2./(1-u).^3};
[x,u,inform] = BVP1DNL(interval,1,0,0,1,f,0,0,u,'Tol',1e-5,'display','off');
figure(3); clf; plot(x,u); xlabel('x'); ylabel('u(x)'); hold on
lambda_ListUp = linspace(lambda0,1.4,20);
N_lambda = 50; d_lambda = diff(lambda_ListUp(1:2))/N_lambda;
u_meshUp = zeros(2*N-1,length(lambda_ListUp));
for jj = 1:length(lambda_ListUp)-1
  lambda = lambda_ListUp(jj)
  u_meshUp(:,jj) = u;
  for kk = 0:N_lambda
    lambda = lambda_ListUp(jj) + kk*d_lambda;
    f = {@(x,u)lambda./(1-u).^2, @(x,u)+lambda*2./(1-u).^3};
    [x,u,inform] = BVP1DNL(interval,1,0,0,1,f,0,0,u,'Tol',1e-5,'display','off');
  endfor % kk
  figure(3); plot(x,u); xlabel('x'); ylabel('u'); drawnow()
endfor
u_meshUp(:,end) = u;
hold off
figure(4); waterfall(x,lambda_ListUp,u_meshUp'); title('upper branch')
           xlabel('x'); ylabel('\lambda'); zlabel('u(x)'); view(-140,35)

endswitch

printing = 0;
if printing
  figure(2); print -dpdfcrop MEMSLowerBranch.pdf
  figure(4); print -dpdfcrop MEMSUpperBranch.pdf
endif

