clear all
pkg load femoctave
E = 110e9; nu = 0.35; %%% copper
%%E = 200e9; nu = 0.25; %%% steel
R = 0.1; dR = 0.01;
nR = 5; nPhi = 51;  %% number of layers in radial and angular direction
Estar = E/(1-nu^2);  nustar = nu/(1-nu);

global P
P = 10e5;  %% 10 atm pressure
FEMmesh = CreateMeshRect(linspace(R,R+dR,nR+1),linspace(0,pi/2,nPhi+1),-21,-12,-33,-22);
function new_xy = Deform(xy)
  new_xy = [xy(:,1).*cos(xy(:,2)),xy(:,1).*sin(xy(:,2))];
endfunction
FEMmesh = MeshDeform(FEMmesh,'Deform');
FEMmesh = MeshUpgrade(FEMmesh,'quadratic');

%% define the radial pressure to be applied on the inside
function res = gN1(xy)
  global P
  angle = atan2(xy(:,2),xy(:,1));
  res   = P*cos(angle);
endfunction
function res = gN2(xy)
  global P
  angle = atan2(xy(:,2),xy(:,1));
  res   = P*sin(angle);
endfunction

[u1,u2] = PlaneStrain(FEMmesh,E,nu,{0,0},{0,0},{'gN1','gN2'});

factor = 400;
figure(1); trimesh(FEMmesh.elem,FEMmesh.nodes(:,1)+factor*u1,FEMmesh.nodes(:,2)+factor*u2,'color','red','linewidth',2);
hold on ;  trimesh(FEMmesh.elem,FEMmesh.nodes(:,1),FEMmesh.nodes(:,2),'color','green','linewidth',1);
hold off; axis equal; xlabel('x'); ylabel('y');

[sigma_x,sigma_y,tau_xy,sigma_z] = EvaluateStress(FEMmesh,u1,u2,E,nu);
vonMises = EvaluateVonMises(sigma_x,sigma_y,tau_xy,sigma_z);

figure(2); FEMtrimesh(FEMmesh,vonMises); xlabel('x'); ylabel('y');
           title('von Mises stress'); view([25,25])
vonMises_min_max = [min(vonMises),max(vonMises)]


%% evaluation at one angle, all radii
alpha = pi/4; Nr = 101; Nmid = (Nr+1)/2; %% use an odd number for Nr
r = linspace(R,R+dR,Nr)'; x = r*cos(alpha); y = r*sin(alpha);

%%[sigma_x,sigma_y,tau_xy,sigma_z] = EvaluateStress(FEMmesh,u1,u2,E,nu);
[sigma_1,sigma_2] = EvaluatePrincipalStress(sigma_x,sigma_y,tau_xy);
sigma_1r = FEMgriddata(FEMmesh,sigma_1,x,y);
sigma_2r = FEMgriddata(FEMmesh,sigma_2,x,y);
sigma_3r = nu*(sigma_1r+sigma_2r);
sigma_2r_min_max = [min(sigma_2r),max(sigma_2r)]

figure(3); plot(r,[sigma_1r,sigma_2r,sigma_3r]);
           xlabel('radius'); ylabel('principal stresses')
           ;legend('\sigma_1','\sigma_2','\sigma_3','location','west')

%% examine stress at middle point
x_mid = x(Nmid); y_mid = y(Nmid);

sigma_x = FEMgriddata(FEMmesh,sigma_x,x_mid,y_mid);
sigma_y = FEMgriddata(FEMmesh,sigma_y,x_mid,y_mid);
tau_xy  = FEMgriddata(FEMmesh,tau_xy ,x_mid,y_mid);
RotMat = [cos(alpha) -sin(alpha);+sin(alpha) cos(alpha)];
stress   = [sigma_x,tau_xy;tau_xy,sigma_y]
stress_rotated   = RotMat'*stress*RotMat

