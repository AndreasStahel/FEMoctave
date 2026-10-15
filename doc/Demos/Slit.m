## -*- texinfo -*-
## @deftypefn  {} {} Slit.m
##
## This is a demo file  inside the `doc/Examples/Elasticity/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

%% construct the domain
Hole.name = "Hole";  W = 0.33;
Hole.border = [0,-W,-22;-5e-3,0,-22;0,+W,-22;+5e-3,0,-22];
Hole.point = [0,0];
Segment.name= "Segment"; Segment.border = [-5e-3,W+5e-3,0;+5e-3,W+5e-3,0];
xy = [-1,-1,-22; 1,-1,-11;+1,+1,-22;-1,+1,-11];
Mesh = CreateMeshTriangle("Slit",xy,2e-2,Hole,Segment);
figure(1); FEMtrimesh(Mesh); axis equal
%Mesh = MeshUpgrade(Mesh,'quadratic');
Mesh = MeshUpgrade(Mesh,'cubic');
figure(2); FEMtrimesh(Mesh); xlim([-0.1, +0.1]); ylim([W-0.1,W+0.1]); axis equal

Displacement = @(xy)0.1*xy(:,1);
E = 1; nu = 0.25;
[u1,u2] = PlaneStress(Mesh,E,nu,{0,0},{Displacement,0},{0,0});  %% solve the plane stress problem

figure(11); FEMtrimesh(Mesh,u1); xlabel('x'); ylabel('y'); zlabel('u_x'); view([-20,10])
figure(12); FEMtrimesh(Mesh,u2); xlabel('x'); ylabel('y'); zlabel('u_y'); view([-20,10])

[eps_xx,eps_yy,eps_xy]   = EvaluateStrain(Mesh,u1,u2);
figure(13); FEMtrimesh(Mesh,eps_xx); xlabel('x'); ylabel('y'); zlabel('\epsilon_{xx}')
figure(14); FEMtrimesh(Mesh,eps_yy); xlabel('x'); ylabel('y'); zlabel('\epsilon_{yy}')

[sigma_x,sigma_y,tau_xy] = EvaluateStress(Mesh,u1,u2,E,nu);
vonMises = EvaluateVonMises(sigma_x,sigma_y,tau_xy);
figure(21); FEMtrisurf(Mesh,vonMises); xlabel('x'); ylabel('y'); zlabel('von Mises')
figure(22); FEMtricontour(Mesh,vonMises,51); xlabel('x'); ylabel('y');
figure(23); FEMtrimesh(Mesh,sigma_x); xlabel('x'); ylabel('y'); zlabel('\sigma_x')
figure(24); FEMtrimesh(Mesh,sigma_y); xlabel('x'); ylabel('y'); zlabel('\sigma_y')
figure(25); FEMtrimesh(Mesh,tau_xy); xlabel('x'); ylabel('y'); zlabel('\tau_{xy}')
Max_vonMises = max(vonMises)


