%% 1 GPa = 10^9 N/m^2 = 10^3 MPa = 10^3 N/mm^2
E = 200e3; nu = 0.24 ;%% N/mm^2
%% parameters for steel E = 200 GPa = 200e3 MPa, nu = 0.24, yield strength = 300 MPa
%% parameters for gold : A 79 GPa = 79e3 MPa, nu = 0.42, yield strength = 200 MPa
R = 1; D = 0.01;    %% radius of cylinder and indentation depth
W = 0.75; H = 1;    %% width and height of the computational domain
FudgeFactor = 0.8469

global a P
Estar = E/(1-nu^2);
P = 4*Estar/3*sqrt(R)*D^(3/2)
a = (3*P*R/(4*Estar))^(1/3)
D = (3*P/(4*Estar))^(2/3)/R^(1/3);
if 0  %% apply the FudgeFactor
  D = D/FudgeFactor
endif

function res = load(rz)
  global a P
   r = rz(:,1);
  res = -3*P/(2*pi*a^2)*sqrt(1-r.^2/a^2).*(r<a);
endfunction

r_low = 1e-7;
N = 51;
if 0  %% mesh by rectangle
  Mesh = CreateMeshRect(linspace(r_low,W,N),linspace(-H,0,N),-11,-23,-12,-12);
else  %% mesh by triangle
  area = 0.5*W*H/N^2;
  dd = 0.001;
  Mesh = CreateMeshTriangle('Flat',[r_low,0,-23;a-dd,0,-23;a,0,-22;a+dd,0,-22;W,0,-12;W,-H,-11;r_low,-H,-12],area);
endif
Mesh = MeshUpgrade(Mesh,'quadratic');
%%Mesh = MeshUpgrade(Mesh,'cubic');
[ur,uz] = AxiStress(Mesh,E,nu,{0,0},{0,0},{0,'load'});

figure(1); FEMtrimesh(Mesh,ur); xlabel('r'); ylabel('z'); zlabel('u_r')
figure(2); FEMtrimesh(Mesh,uz); xlabel('r'); ylabel('z'); zlabel('u_z')
r = linspace(r_low,2*a,1000)';
uz_edge = FEMgriddata(Mesh,uz,r,zeros(size(r)));
circle = R-D-sqrt(R^2-r.^2);
figure(11); plot(r,uz_edge,'r',r,circle,'k',[a,a],[-D,0],'g')
            xlabel('r'); ylabel('u_z'); xlim([0,2*a]);
            legend('u_z','circle','r=a','location','northwest')
DifferenceUzAtOrigin = circle(1)-uz_edge(1)

[sigma_r,sigma_y,sigma_z,tau_xz] = EvaluateStressAxi(Mesh,ur,uz,E,nu);
VonMises = EvaluateVonMisesAxi(sigma_r,sigma_y,sigma_z,tau_xz);

NN = 51;
[rr,zz] = meshgrid(linspace(r_low,3*a,NN),linspace(-3*a,0,NN));
sigma_zg = FEMgriddata(Mesh,sigma_z,rr,zz);
sigma_rg = FEMgriddata(Mesh,sigma_r,rr,zz);
VonMises_g = FEMgriddata(Mesh,VonMises,rr,zz);
figure(31); contourf(rr,zz,sigma_rg/1e3,51); xlabel('r'); ylabel('z');
            title('\sigma_r [kPa]'); axis equal; colorbar
figure(32); contourf(rr,zz,sigma_zg/1e3,51); xlabel('r'); ylabel('z');
            title('\sigma_z [kPa]'); axis equal; colorbar
figure(33); contourf(rr,zz,VonMises_g/1e3,51); xlabel('r'); ylabel('z');
            title('von Mises [kPa]'); axis equal; colorbar

