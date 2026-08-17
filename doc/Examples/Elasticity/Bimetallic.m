## -*- texinfo -*-
## @deftypefn  {} {} Bimetallic.m
##
## This is a demo file  inside the `doc/Examples/Elasticity/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

%% horizontal split of materials
E1 = 120e3; nu1 = 0.34; alpha1 = 16e-6; %% copper
E2 = 210e3; nu2 = 0.30; alpha2 = 11e-6; %% steel
L = 40; D1 = 0.6; D2 = 0.4; DT = 70;
%%nu1 = 0; nu2 = 0;
E  = @(xy) E1*(xy(:,2)<0) +  E2*(xy(:,2)>= 0);
nu = @(xy)nu1*(xy(:,2)<0) + nu2*(xy(:,2)>= 0);
ThermalCoeff = @(xy)DT*alpha1*(xy(:,2)<0) + DT*alpha2*(xy(:,2)>= 0);

Mesh = CreateMeshRect(linspace(0,L,41),linspace(-D1,D2,41),-22,-22,-12,-22);
Mesh = MeshAddConstraint(Mesh,[0,0],[-1,-1]);
Mesh = MeshUpgrade(Mesh,'cubic');

[u1,u2] = PlaneStress(Mesh,E,nu,{0,0},{0,0},{0,0},'thermal',ThermalCoeff);
figure(1); FEMtrimesh(Mesh,u1); xlabel('x'); ylabel('y'); zlabel('u_1')
figure(2); FEMtrimesh(Mesh,u2); xlabel('x'); ylabel('y'); zlabel('u_2')
Max_u1u2 = [max(u1),max(u2)]

[sigma_x,sigma_y,tau_xy] = EvaluateStress(Mesh,u1,u2,E,nu,'thermal',ThermalCoeff);
figure(3); FEMtrimesh(Mesh,sigma_x); xlabel('x'); ylabel('y'); zlabel('\sigma_x');
           zlim([-30,50]); view([-125,25])

y = linspace(-D1,D2,1000)'; xy = [L/2*ones(size(y)),y];
[sigma_xc,sigma_yc,tau_xyc] = ...
EvaluateStress(Mesh,u1,u2,E,nu,'curve',xy,'thermal',ThermalCoeff);
figure(4); plot(y,sigma_xc); xlabel('y'); ylabel('\sigma_x')

x = linspace(0,L)';
u2_x = FEMgriddata(Mesh,u2,x,zeros(size(x)));
figure(5); plot(x,u2_x); xlabel('x'); ylabel('u_2(x,0)')
Param_u2 = polyfit(x,u2_x,2)
Curvature_FEM = 2*Param_u2(1)
