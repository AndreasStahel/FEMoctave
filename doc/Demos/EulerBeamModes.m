## -*- texinfo -*-
## @deftypefn  {} {} EulerBeamModes.m
##
## This is a demo file  inside the `doc/Demos/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

clear *
L = 0.20; H = 0.01; W = 0.01; rho = 2.7e3;
E = 70e9; nu = 0.33; %% Aluminum
I2 = 1/12*H^3*W;
Mode = 2
Nx = 20; Ny = 3;
MeshType = 'linear';
MeshType = 'quadratic';
%MeshType = 'cubic';

f = @(z) 1+cos(z).*cosh(z);  %% clamped at x=0, free at x=L
% z = linspace(0,Mode*pi,100);
% figure(101); plot(z,f(z)); xlabel('z'); ylabel('f(z)')
z0 = fsolve(f,Mode*pi-pi/2);
freqEuler = z0^2*sqrt(E*I2/(rho*H*W))/(2*pi*L^2)
% p = k*lambda^0.25*L;
% C = (sin(p)-sinh(p))/(cos(p)+cosh(p));
% x = linspace(0,L);  x_p = k*lambda^0.25*x;
% y = cos(x_p)-cosh(x_p) + C*(sin(x_p)-sinh(x_p)); y = y/y(end);
% figure(102); plot(x,y); xlabel('x'); ylabel('height y(x)');
Mesh = CreateMeshRect(linspace(0,L,Nx+1),linspace(0,+H,Ny+1),-22,-22,-11,-22);
switch MeshType
 case 'quadratic'
   Mesh = MeshUpgrade(Mesh,'quadratic');
 case 'cubic'
   Mesh = MeshUpgrade(Mesh,'cubic');
endswitch
[la,u1,u2] = PlaneStressEig(Mesh,E,nu,rho,max(4,Mode));
freqFEM = sqrt(la')/(2*pi)
u1_disp = u1(:,Mode); u2_disp = u2(:,Mode);
MaxDisp = max(max(abs(u1_disp)),max(abs(u2_disp)));
u1_disp = u1_disp/MaxDisp; u2_disp = u2_disp/MaxDisp;
figure(1);FEMtrimesh(Mesh,u1_disp); xlabel('x'); ylabel('y'); zlabel('u_1')
figure(2);FEMtrimesh(Mesh,u2_disp); xlabel('x'); ylabel('y'); zlabel('u_2')
[sigma_x,sigma_y,tau_xy] = EvaluateStress(Mesh,u1_disp,u2_disp,E,nu);
figure(11); FEMtrimesh(Mesh,sigma_x*1e-9); xlabel('x'); ylabel('y'); zlabel('\sigma_x [GPa]')
figure(12); FEMtrimesh(Mesh,sigma_y*1e-9); xlabel('x'); ylabel('y'); zlabel('\sigma_y [GPa]')
figure(13); FEMtrimesh(Mesh,tau_xy*1e-9);  xlabel('x'); ylabel('y'); zlabel('\tau_{xy} [GPa]')

figure(20);clf
factor = 0.005;
trimesh(Mesh.elem,Mesh.nodes(:,1)+factor*u1_disp,Mesh.nodes(:,2)+factor*u2_disp,...
'color','red','linewidth',2);
hold on ;
trimesh(Mesh.elem,Mesh.nodes(:,1),Mesh.nodes(:,2),'color','green','linewidth',1);
xlabel('x'); ylabel('y'); %%xlim([0,L*1.1])

