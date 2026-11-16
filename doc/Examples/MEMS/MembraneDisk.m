CASE = 1;

switch CASE
  case 1  %% construct the lower branch
R = 0.5;
N = 51; interval = linspace(0,R,N)';
lambda_List = [0.0:0.1:3.15];
lambda = lambda_List(1);
a = @(r) r;
f = {@(r,u)lambda*r./(1-u).^2 , @(r,u)+lambda*2*r./(1-u).^3};
[r,u,inform] = BVP1DNL(interval,a,0,0,1,f,[0,0],0,0,'Tol',1e-5,'display','iter');
figure(1); clf; plot(r,u); xlabel('r'); ylabel('u')
hold on
u_mesh = zeros(2*N-1,length(lambda_List));
u_mesh(:,1) = u;
for jj = 2:length(lambda_List)
  lambda = lambda_List(jj)
  f = {@(r,u)lambda*r./(1-u).^2, @(r,u)+lambda*2*r./(1-u).^3};
  [r,u,inform] = BVP1DNL(interval,a,0,0,1,f,[0,0],0,u,'Tol',1e-5,'display','iter');
  figure(1); plot(r,u); xlabel('r'); ylabel('u')
  u_mesh(:,jj) = u;  u0_List(jj) = u(1);
endfor
hold off
figure(2); waterfall(r,lambda_List,u_mesh'); title('lower branch');
           xlabel('r'); ylabel('\lambda'); zlabel('u(r)'); view(40,45); ylim('tight')
figure(11); plot(lambda_List,u_mesh(1,:)); xlabel('\lambda'); ylabel('u(0)')

  case 2  %% construct the upper branch
load UpperSolutionDisk.mat;  r_max = r; u_max = u;
R = 0.5; N = 51; interval = linspace(0,R,N)';
r = sort([interval; interval(1:end-1) + diff(interval)/2]);
u = interp1(r_max,u_max,r); u0 = u; u0(1) = u0(2); lambda0 = lambda;
a = @(r) r;
f = {@(r,u)lambda*r./(1-u).^2, @(r,u)+lambda*2*r./(1-u).^3};
[r,u,inform] = BVP1DNL(interval,a,0,0,1,f,[0,0],0,u0,'Tol',1e-5,'display','off');
figure(3); clf; plot(r,u); xlabel('r'); ylabel('u(r)'); hold on
lambda_ListUp = linspace(lambda0,3.15,20);
N_lambda = 50; d_lambda = diff(lambda_ListUp(1:2))/N_lambda;
u_meshUp = zeros(2*N-1,length(lambda_ListUp));
for jj = 1:length(lambda_ListUp)-1
  lambda = lambda_ListUp(jj)
  u_meshUp(:,jj) = u;
  for kk = 0:N_lambda
    lambda = lambda_ListUp(jj) + kk*d_lambda;
    f = {@(r,u)lambda*r./(1-u).^2, @(r,u)+lambda*2*r./(1-u).^3};
    [r,u,inform] = BVP1DNL(interval,a,0,0,1,f,[0,0],0,u,'Tol',1e-5,'display','off');
  endfor % kk
  figure(3); plot(r,u); xlabel('r'); ylabel('u'); drawnow()
endfor
u_meshUp(:,end) = u;
hold off
figure(4); waterfall(r,lambda_ListUp,u_meshUp'); title('upper branch')
           xlabel('r'); ylabel('\lambda'); zlabel('u(r)'); view(45,35); ylim('tight')
endswitch

