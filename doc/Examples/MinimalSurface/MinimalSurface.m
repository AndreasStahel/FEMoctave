## -*- texinfo -*-
## @deftypefn  {} {} MinimalSurface.m
##
## This is a demo file  inside the `doc/Examples/MinimalSurface/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

clear *
%%FEMmesh = CreateMeshRect(linspace(-1,1,21),linspace(-1,1,21),-1, -1, -1, -1);
xy = [1,0,-1;0,1,-1;-1,0,-1;0,-1,-1];
FEMmesh = CreateMeshTriangle("square",xy,0.01);
%FEMmesh = MeshUpgrade(FEMmesh);

x = FEMmesh.nodes(:,1); y = FEMmesh.nodes(:,2); elem = FEMmesh.elem;
function res = BC(xy)
  res = abs(xy(:,1));
endfunction

u = BVP2Dsym(FEMmesh,1,0,0,'BC',0,0);
difference = zeros(5,1); area = difference;
for ii = 1:5
  [~,grad] = FEMEvaluateGP(FEMmesh,u);
  coeff = sqrt(1+grad(:,1).^2+ grad(:,2).^2);
  area(ii) = FEMIntegrate(FEMmesh,coeff);
  u_new = BVP2Dsym(FEMmesh,coeff,0,0,'BC',0,0);
  difference(ii) = mean(abs(u_new-u));
  u = u_new;
endfor

Area_Difference = [area,difference]

figure(1)
FEMtrisurf(FEMmesh,u)
xlabel('x'); ylabel('y'); zlabel('z')
