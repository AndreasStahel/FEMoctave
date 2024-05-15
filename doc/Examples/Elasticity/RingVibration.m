## -*- texinfo -*-
## @deftypefn  {} {} RingVibration.m
##
## This is a demo file  inside the `doc/Examples/Elasticity/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

clear *
R = 0.03; D = 0.002;
Nring = 36; angles = linspace(0,2*pi,Nring+1)'; angles = angles(1:end-1);

Ring = [(R+D/2)*cos(angles),(R+D/2)*sin(angles),-22*ones(size(angles))];
Hole.name   = 'Hole';
Hole.border = [(R-D/2)*cos(angles),(R-D/2)*sin(angles),-22*ones(size(angles))];
Hole.point  = [0,0.01];

FEMmesh = CreateMeshTriangle('Ring',Ring,5e-7, Hole);
FEMmesh = MeshUpgrade(FEMmesh,'quadratic');

E = 200e9; nu = 0.25; rho = 8e3; %% steel
Nmodes = 12;
[lambda,u1_all,u2_all] = PlaneStressEig(FEMmesh,E,nu,rho,Nmodes);
frequencies = sqrt(abs(lambda'))/(2*pi)

for Mode = 1:Nmodes
  u1 = u1_all(:,Mode); u2 = u2_all(:,Mode); scale = 3e-3/(max(abs([u1;u2])));
  figure(10+Mode); ShowDeformation(FEMmesh,u1,u2,scale); axis equal; title(sprintf('mode %i',Mode))
endfor



