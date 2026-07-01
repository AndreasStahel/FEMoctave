pkg load femoctave
%% set the parameters
E = 110e3; nu = 0.35; alpha = 1.7e-5; DeltaT = 6; %% copper, length in mm
L = 1;  N = 11; R = L/2;
ThermalCoeff = @(xy)(DeltaT*alpha*( (xy(:,1).^2 + xy(:,2).^2) < R^2));

%% generate the mesh with high resolution along the circle
phi = linspace(0,2*pi,360/3)'; phi = phi(2:end);
Seg.name   = 'Segment';
Seg.border = [R*cos(phi),R*sin(phi),zeros(size(phi))];
BC = -22;
Mesh = CreateMeshTriangle('block',[-L,-L,BC;+L,-L,BC;+L,+L,BC;-L,+L,BC],(L/N)^2,Seg);
Mesh = MeshAddConstraint(Mesh,[0,0],[-1,-1]); Mesh = MeshAddConstraint(Mesh,[R,0],[-2,-1]);
figure(1); FEMtrimesh(Mesh); axis equal
% Mesh = MeshUpgrade(Mesh,'quadratic');
Mesh = MeshUpgrade(Mesh,'cubic');

%% solve the problem
[u1,u2] = PlaneStress(Mesh,E,nu,{0,0},{0,0},{0,0},'thermal',ThermalCoeff);
figure(2); FEMtrimesh(Mesh,u1); xlabel('x'); ylabel('y'); zlabel('u_1')
figure(3); FEMtrimesh(Mesh,u2); xlabel('x'); ylabel('y'); zlabel('u_2')

%% display strains and stresses
[eps_xx,eps_yy,eps_xy] = EvaluateStrain(Mesh,u1,u2);
figure(4); FEMtrimesh(Mesh,eps_xx); xlabel('x'); ylabel('y'); zlabel('\epsilon_{xx}')
figure(5); FEMtrimesh(Mesh,eps_yy); xlabel('x'); ylabel('y'); zlabel('\epsilon_{yy}')

[sigma_x,sigma_y,tau_xy] = EvaluateStress(Mesh,u1,u2,E,nu,'thermal',ThermalCoeff);
figure(6); FEMtrimesh(Mesh,sigma_x); xlabel('x'); ylabel('y'); zlabel('\sigma_x')
figure(7); FEMtrimesh(Mesh,sigma_y); xlabel('x'); ylabel('y'); zlabel('\sigma_y')
figure(8); FEMtrimesh(Mesh,tau_xy);  xlabel('x'); ylabel('y'); zlabel('\tau_{xy}')

%% evaluate along the circle
phi = linspace(0,2*pi,2001)'; phiD = phi/pi*180; Rc = R*1.001;
x = Rc*cos(phi); y = Rc*sin(phi);
[sigma_xc,sigma_yc,tau_xyc] = EvaluateStress(Mesh,u1,u2,E,nu,'thermal',ThermalCoeff,'curve',[x,y]);
sigma_nc = sigma_xc.*cos(phi).^2 + 2*tau_xyc.*cos(phi).*sin(phi) + sigma_yc.*sin(phi).^2;
figure(9); plot(phiD,sigma_xc,phiD,sigma_yc,phiD,sigma_nc); xlabel('angle'); ylabel('\sigma')
           legend('\sigma_x','\sigma_y','\sigma_n','location','northeast'); xlim([0,360])

printing = 0
if printing
  figure(1); print -dpng CopperThermalMesh.png
  figure(2); print -dpng CopperThermal_u1.png
  figure(3); print -dpng CopperThermal_u2.png
  figure(4); print -dpng CopperThermal_epsxx.png
  figure(5); print -dpng CopperThermal_epsyy.png
  figure(9); print -dpdfcrop CopperThermal_sigma.pdf
endif

