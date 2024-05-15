## -*- texinfo -*-
## @deftypefn  {} {} WaiitMitchell.m
##
## This is a demo file  inside the `doc/Examples/Elasticity/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

E = 3e6; nu = 0.45;    %% rubber
Dom = [0,1/3,-11;1/3,1/3,-11;1/3,0,-22;1,0,-22;1,2/3,-11;2/3,1,-22;0,1,-22];
Mesh = CreateMeshTriangle('Domain',Dom,1e-3);
figure(1); FEMtrimesh(Mesh); axis equal
Mesh = MeshUpgrade(Mesh,'cubic');

function res = Disp(xy)
  res = 0.3 .*(sum(xy,2)>1) ; %% if x+y>1
endfunction

[u1,u2] = PlaneStress(Mesh,E,nu,{0,0},{'Disp','Disp'},{0,0});
figure(2); FEMtrisurf(Mesh,u1); xlabel('x'); ylabel('y'); zlabel('u_1')
figure(3); FEMtrisurf(Mesh,u2); xlabel('x'); ylabel('y'); zlabel('u_2')
figure(100); ShowDeformation(Mesh,u1,u2,1); axis equal

[sigma_x,sigma_y,tau_xy] = EvaluateStress(Mesh,u1,u2,E,nu);
figure(4); FEMtrimesh(Mesh,sigma_x); xlabel('x'); ylabel('y'); zlabel('\sigma_x')
figure(5); FEMtrimesh(Mesh,sigma_y); xlabel('x'); ylabel('y'); zlabel('\sigma_y')
figure(6); FEMtrimesh(Mesh,tau_xy); xlabel('x'); ylabel('y'); zlabel('\tau_{xy}')

VonMises = EvaluateVonMises(sigma_x,sigma_y,tau_xy);
[eps_xx,eps_yy,eps_xy] = EvaluateStrain(Mesh,u1,u2);
EnergyDensity = EvaluateEnergyDensity(Mesh,eps_xx,eps_yy,eps_xy,E,nu);

ElasticEnergy = FEMIntegrate(Mesh,EnergyDensity)
Force = 2*ElasticEnergy/(sqrt(2)*0.3)
Contour = [Dom(:,1:2);Dom(1,1:2)];
figure(7); clf; FEMtricontour(Mesh,VonMises/1e6,[0:0.1:2]); xlabel('x'); ylabel('y');
           hold on; plot(Contour(:,1),Contour(:,2),'k'); colorbar();
figure(8); clf; FEMtricontour(Mesh,EnergyDensity/1e6,[0:0.05:1]); xlabel('x'); ylabel('y');
           hold on; plot(Contour(:,1),Contour(:,2),'k'); colorbar();

eps_45 = (eps_xx+eps_yy+2*eps_xy)/2;
x = linspace(1/3,1-1/6)';
eps_45line = FEMgriddata(Mesh,eps_45,x,x);
figure(9); clf; FEMtricontour(Mesh,eps_45,[0:0.05:0.6]); xlabel('x'); ylabel('y');
           hold on; plot(Contour(:,1),Contour(:,2),'k'); colorbar();
figure(10); plot(x,eps_45line); xlabel('x'); ylabel('eps_{45}')
Stretch = [trapz(sqrt(2)*x,eps_45line),sqrt(2)*0.3]
