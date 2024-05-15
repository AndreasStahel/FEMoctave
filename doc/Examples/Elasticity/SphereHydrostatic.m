## -*- texinfo -*-
## @deftypefn  {} {} SphereHydrostatic.m
##
## This is a demo file  inside the `doc/Examples/Elasticity/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn


clear *
global P
P = -0.1;
E = 1; nu = 0.25; R = 1;
N = 51; alpha = linspace(0,+pi/2,N)'; x = R*cos(alpha); z = R*sin(alpha);
Dom = [0,0,-21;[x,z,-33*ones(size(x))]]; Dom(end,3) = -12;
function res = gNr(rz)
  global P
  alpha = atan2(rz(:,2),rz(:,1));
  res   = P*cos(alpha);
endfunction
function res = gNz(rz)
  global P
  alpha = atan2(rz(:,2),rz(:,1));
  res   = P*sin(alpha);
endfunction

Mesh = CreateMeshTriangle('quart',Dom,1e-2);
Mesh = MeshUpgrade(Mesh,'quadratic');
[ur,uz] = AxiStress(Mesh,E,nu,{0,0},{0,0},{'gNr','gNz'});
figure(1); FEMtrimesh(Mesh,ur); xlabel('r'); ylabel('z'); zlabel('u_r'); view([20,40])
figure(2); FEMtrimesh(Mesh,uz); xlabel('r'); ylabel('z'); zlabel('u_z'); view([-120,20])

[eps_xx,eps_yy,eps_zz,eps_xz] = EvaluateStrainAxi(Mesh,ur,uz);
figure(3); FEMtrimesh(Mesh,eps_xx); xlabel('r'); ylabel('z'); zlabel('eps_{xx}')
figure(4); FEMtrimesh(Mesh,eps_yy); xlabel('r'); ylabel('z'); zlabel('eps_{yy}')
figure(5); FEMtrimesh(Mesh,eps_zz); xlabel('r'); ylabel('z'); zlabel('eps_{zz}')
figure(6); FEMtrimesh(Mesh,eps_xz); xlabel('r'); ylabel('z'); zlabel('eps_{xz}')
W = EvaluateEnergyDensityAxi(Mesh,eps_xx,eps_yy,eps_zz,eps_xz,E,nu);
figure(7); FEMtrimesh(Mesh,W); xlabel('r'); ylabel('z'); zlabel('energy density')
r = Mesh.nodes(:,1);
EnergyIntegrated = FEMIntegrate(Mesh,2*pi*r.*W)
EnergyDensity    = EnergyIntegrated/(4/3*pi*R^3/2)
