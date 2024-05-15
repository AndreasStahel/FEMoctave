## -*- texinfo -*-
## @deftypefn  {} {} PlaneStrainExample ()
##
## This is a demo file  inside the `doc/Demos/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

W = 0.1; H = 0.3; E = 110e9; nu = 0.35;  %% copper

FEMmesh = CreateMeshRect(linspace(-W/2,W/2,10),linspace(0,H,30),-11,-12,-22,-22);
function xy_new = Deform(xy)
  xy_new = [xy(:,1).*(1-0.5/0.3*xy(:,2)) , xy(:,2)];
endfunction
FEMmesh = MeshDeform(FEMmesh,'Deform');
CMesh = CreateMeshRect(linspace(-W/2,W/2,6),linspace(0,H,20),-11,-12,-22,-22);
CMesh = MeshDeform(CMesh,'Deform'); %% create a course mesh on the same domain

figure(1); FEMtrimesh(FEMmesh)
           xlabel('x'); ylabel('y'); axis equal

%%FEMmesh = MeshUpgrade(FEMmesh,'quadratic');  %% uncomment for second order elements
FEMmesh = MeshUpgrade(FEMmesh,'cubic');  %% uncomment for third order elements

f = {0,0}; gN = {0,0};
function res = gD(xy)
  res = +(xy(:,2)>0.1)*0.01;
endfunction

[u1,u2] = PlaneStrain(FEMmesh,E,nu,f,{'gD',0},gN);

u1i = FEMgriddata(FEMmesh,u1,CMesh.nodes(:,1),CMesh.nodes(:,2));
u2i = FEMgriddata(FEMmesh,u2,CMesh.nodes(:,1),CMesh.nodes(:,2));
figure(1); ShowDeformation(CMesh,u1i,u2i,2);
           axis equal; xlabel('x'); ylabel('y'); ylim([0,0.35])
figure(2); FEMtrimesh(FEMmesh,u1)
           xlabel('x'); ylabel('y'); zlabel('u_1'); view([50,30])
figure(3); FEMtrimesh(FEMmesh,u2)
           xlabel('x'); ylabel('y'); zlabel('u_2'); view([50,30])


[eps_xx,eps_yy,eps_xy] = EvaluateStrain(FEMmesh,u1,u2);
figure(4);
subplot(1,3,1); FEMtrimesh(FEMmesh,eps_xx)
                xlabel('x'); ylabel('y'); zlabel('\epsilon_{xx}'); view([50,30])
subplot(1,3,2); FEMtrimesh(FEMmesh,eps_yy)
                xlabel('x'); ylabel('y'); zlabel('\epsilon_{yy}'); view([50,30])
subplot(1,3,3); FEMtrimesh(FEMmesh,eps_xy)
                xlabel('x'); ylabel('y'); zlabel('\epsilon_{xy}'); view([50,30])

[sigma_x,sigma_y,tau_xy,sigma_z] = EvaluateStress(FEMmesh,u1,u2,E,nu);

figure(5); title('stress')
subplot(1,3,1); FEMtrimesh(FEMmesh,sigma_x)
                xlabel('x'); ylabel('y'); zlabel('\sigma_x'); view([50,30])
subplot(1,3,2); FEMtrimesh(FEMmesh,sigma_y)
                xlabel('x'); ylabel('y'); zlabel('\sigma_y'); view([50,30])
subplot(1,3,3); FEMtrimesh(FEMmesh,tau_xy)
                xlabel('x'); ylabel('y'); zlabel('\tau_{xy}'); view([50,30])

vonMises = EvaluateVonMises(sigma_x,sigma_y,tau_xy,sigma_z);
figure(6); FEMtrimesh(FEMmesh,vonMises)
           xlabel('x'); ylabel('y'); zlabel("von Mises stress"); view([120,30])

Tresca = EvaluateTresca(sigma_x,sigma_y,tau_xy,sigma_z);
figure(7); FEMtrimesh(FEMmesh,Tresca)
           xlabel('x'); ylabel('y'); zlabel('Tresca stress'); view([120,30])

[s1,s2] = EvaluatePrincipalStress(sigma_x,sigma_y,tau_xy);
