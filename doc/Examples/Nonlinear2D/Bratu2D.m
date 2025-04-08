## -*- texinfo -*-
## @deftypefn  {} {} Bratu2D.m
##
## This is a example file inside the `doc/Examples/Nonlinear2D/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

CASE = 1;
MeshCase  = 2;  N = 10;
switch MeshCase
case 1  %% full square
  Mesh = CreateMeshRect(linspace(0,1,2*N+1),linspace(0,1,2*N+1),-1,-1,-1,-1);
case 2  %% quarter of the square
  Mesh = CreateMeshRect(linspace(0,0.5,N+1),linspace(0,0.5,N+1),-1,-2,-1,-2);
case 3  %% quarter of the square, triangle
  Mesh = CreateMeshTriangle('Shape',[0,0,-1;0.5,0,-2;0.5,0.5,-2;0,0.5,-1],0.2/N^2);
case 4  %% eigth of the square, triangle
  Mesh = CreateMeshTriangle('Shape',[0,0,-2;0.5,0.5,-2;0,0.5,-1],0.2/N^2);
endswitch
Mesh = MeshUpgrade(Mesh,'cubic');

switch CASE
case 1  %% BVP2DNL, one solution, lower branch
  C = 1.0;
  u0 = @(xy)0.8*sin(pi*xy(:,1)).*sin(pi*xy(:,2));
  f = {@(x,u)C*exp(u),@(x,u)C*exp(u)};
  u = BVP2DNL(Mesh,1,0,0,0,f,0,0,0,u0,'display','iter');
  figure(1); FEMtrimesh(Mesh,u)
  u_max = [max(u),FEMgriddata(Mesh,u,0.5,0.5)]

case 2  %% BVP2DNL, one solution, upper branch
  C = 1.0;
  u0 = @(xy)5.5*sin(pi*xy(:,1)).*sin(pi*xy(:,2));
  f = {@(x,u)C*exp(u),@(x,u)C*exp(u)};
  u = BVP2DNL(Mesh,1,0,0,0,f,0,0,0,u0,'maxiter',20);
  figure(2); FEMtrimesh(Mesh,u)
  u_max = [max(u),FEMgriddata(Mesh,u,0.5,0.5)]

case 3  %% BVP2DNL, parameter moving up, lower branch
  C = 1.0;
  u0 = @(xy)0.8*sin(pi*xy(:,1)).*sin(pi*xy(:,2));
  f = {@(x,u)C*exp(u),@(x,u)C*exp(u)};
  u = BVP2DNL(Mesh,1,0,0,0,f,0,0,0,u0,'maxiter',20);
  C_list = C;  u_max_list = max(u);
  step = 0.2; half_counter = 0; iterations = 0;
  do
    u0 = u; C = C+step;
    f = {@(x,u)C*exp(u),@(x,u)C*exp(u)};
    [u,inform] = BVP2DNL(Mesh,1,0,0,0,f,0,0,0,u0,'maxiter',20);
    iterations++;
    if (inform.info == 1)&&(inform.iter<5)
      u_max_list = [u_max_list,max(u)]; C_list = [C_list,C];
      disp(sprintf("%i, C = %f, max(u) = %f, iterations = %i",...
                   iterations,C,u_max_list(end),inform.iter))
      u0 = u;
    else
      C = C-step;
      step = step/2; half_counter++;
    endif
  until or(iterations>=150,half_counter>15)
  figure(3); plot(C_list,u_max_list,'.-');
             xlabel('parameter C'); ylabel('u_\infty')
  C_list3 = C_list; u_max_list3 = u_max_list;
case 4  %% BVP2DNL, parameter moving up, upper branch
  C = 1.0;
  u0 = @(xy)5.5*sin(pi*xy(:,1)).*sin(pi*xy(:,2));
  f = {@(x,u)C*exp(u),@(x,u)C*exp(u)};
  [u,inform] = BVP2DNL(Mesh,1,0,0,0,f,0,0,0,u0,'maxiter',20);
  C_list = C;  u_max_list = max(u);
  disp(sprintf("%i, C=%f, max(u)= %f, iterations = %i",...
                0,C,u_max_list(end),inform.iter))
  step = 0.025; half_counter = 0; iterations = 0;
  do
    u0 = u; C = C+step;
    f = {@(x,u)C*exp(u),@(x,u)C*exp(u)};
    [u,inform] = BVP2DNL(Mesh,1,0,0,0,f,0,0,0,u0,'maxiter',20);
    iterations++;
    if (inform.info == 1)&&(inform.iter<5)
      u_max_list = [u_max_list,max(u)]; C_list = [C_list,C];
      disp(sprintf("%i, C = %f, max(u) = %f, iterations = %i",...
                    iterations,C,u_max_list(end),inform.iter))
      u0 = u;
    else
      C = C-step;
      step = step/2; half_counter++;
    endif
  until or((iterations>=600),(half_counter>15))
  figure(4); plot(C_list,u_max_list,'.-');
             xlabel('parameter C'); ylabel('u_\infty')
  C_list4 = C_list; u_max_list4 = u_max_list;

case 5
  figure(5); plot(C_list3,u_max_list3,C_list4,u_max_list4);
             xlabel('parameter C'); ylabel('u_\infty'); xlim([0,7])
endswitch






