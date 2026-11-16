
CASE = 1 ; %% 1: basic solution, 3: with the geometric nonlinearity

switch CASE
  case 1
L = 1.0; H = 0.1; N1 = 21; N2 = 3;
Mesh = CreateMeshRect(linspace(0,L,N1+1),linspace(-H,0,N2+1),-22,-22,-11,-22);
Mesh = MeshUpgrade(Mesh,'cubic');

lambda = 1.56e3;
fx = 0; fy = +lambda; E = 1e+6; nu = 0.3;
u2atL = 3/2*lambda*L^4/(E*H^2)  %% the elementary Euler beam result

[u1,u2] = PlaneStress(Mesh,E,nu,{fx,fy},{0,0},{0,0});
Max_u2 = max(u2)

for jj = 1:20
  y_GP = Mesh.GP(:,2) + FEMEvaluateGP(Mesh,u2);
  fy = lambda./(1-y_GP).^2;
  [u1,u2] = PlaneStress(Mesh,E,nu,{fx,fy},{0,0},{0,0});
  Max_u2 = max(u2)
endfor
figure(1); FEMtrimesh(Mesh,u2); xlabel('x'); ylabel('y'); zlabel('u_2')
figure(2); ShowDeformation(Mesh,u1,u2,1); axis equal

  case 2
L = 1; H = 0.1; N1 = 21; N2 = 3;
Mesh = CreateMeshRect(linspace(0,L,N1+1),linspace(-H,0,N2+1),-22,-22,-11,-22);
Mesh = MeshUpgrade(Mesh,'cubic');

lambda = 1.56e3;
fx = 0; fy = +lambda; E = 1e+6; nu = 0.3;

[u1,u2] = PlaneStress(Mesh,E,nu,{fx,fy},{0,0},{0,0});
Max_u2 = max(u2)
for jj = 1:10
  [u2_GP, du2_GP] = FEMEvaluateGP(Mesh,u2); Angle = atan(du2_GP(:,1));
  y_GP = Mesh.GP(:,2) + FEMEvaluateGP(Mesh,u2);
  fy = lambda./(1-y_GP).^2.*cos(Angle); fx = lambda./(1-y_GP).^2.*sin(Angle);
  [u1,u2] = PlaneStress(Mesh,E,nu,{fx,fy},{0,0},{0,0});
  Max_u2 = max(u2)
endfor
figure(3); ShowDeformation(Mesh,u1,u2,1);
           xlabel('x'); ylabel('y'); axis equal
endswitch
