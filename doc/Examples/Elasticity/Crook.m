## -*- texinfo -*-
## @deftypefn  {} {} Crook.m
##
## This is a demo file  inside the `doc/Examples/Elasticity/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

pkg load femoctave
W = 0.01; H = 0.1; Load = 1e6;;
Layers = 2*5; gap = W/5;

if 0 %% no rounding
  Domain = [-W -W -22; -W H -11; 0 H -22; 0 gap -22; 0 0 -22; gap 0 -22; H 0 -23; H -W -22];
else %% with a rounded corner
  Domain = [-W -W -22; -W H -11; 0 H -22; 0 gap -22; gap*0.366 gap*0.366 -22; gap 0 -22; H 0 -23; H -W -22];
endif

FEMmesh = CreateMeshTriangle ('Crook1',Domain,(W/Layers)^2);
figure(1); FEMtrimesh(FEMmesh); xlabel('x'); ylabel('y'); axis([-W 3*gap -W 3*gap])
if 1
  FEMmesh = MeshUpgrade(FEMmesh,'quadratic');
else
  FEMmesh = MeshUpgrade(FEMmesh,'cubic');
endif

E = 200e9; nu = 0.25; %%% steel
[u1,u2] = PlaneStress(FEMmesh,E,nu,{0,0},{0,0},{0,-Load});

MaximalDisplacement = min(u2)

[~,slope_x,~] = FEMgriddata(FEMmesh,u2,0,-W/2)
yi = linspace(-0.01,0.1)'; xi =-0.005*ones(size(yi));
u1i = FEMgriddata(FEMmesh,u1,xi,yi);
figure(8); plot(yi,u1i)
           xlabel('y'); ylabel('u_1')

p = LinearRegression([yi.^2,yi,ones(size(yi))],u1i);  %% linear regression of a polynomial of degree 2
slope = polyval([2*p(1) p(2)],-W/2)                   %% evaluate the derivative of the polynomial

CoarseMesh = CreateMeshRect([-W:W/3:H],[-W:W/3:H],-11,-11,-11,-11);
x = CoarseMesh.nodes(:,1); y = CoarseMesh.nodes(:,2);
u1i = FEMgriddata(FEMmesh,u1,x,y); u2i = FEMgriddata(FEMmesh,u2,x,y);
x(isnan(u1i)) = NaN;

factor = H/10/abs(min(u2));
figure(2); clf;
  trimesh(CoarseMesh.elem,x,y,'color','green','linewidth', 1)
hold on;     trimesh(CoarseMesh.elem,x+factor*u1i,y+factor*u2i,'color','red','linewidth', 1)
axis equal; hold off


[sigma_x,sigma_y,tau_xy] = EvaluateStress(FEMmesh,u1,u2,E,nu);

dist = linspace(-W,0,100)'; HH = H/2*ones(size(dist));
sigma_y_slice_H = FEMgriddata(FEMmesh,sigma_y,dist,HH);
figure(3); plot(dist,sigma_y_slice_H/1e6);
           xlabel('x'); ylabel('\sigma_y [MPa]');xlim([-W,0])
Integral_sigma_y = W*trapz(dist,sigma_y_slice_H)
Integral_Moment = W*trapz(dist,dist.*sigma_y_slice_H)

sigma_x_slice_V = FEMgriddata(FEMmesh,sigma_x,HH,dist);
tau_xy_slice_V  = FEMgriddata(FEMmesh,tau_xy,HH,dist);
figure(4); plot(dist,sigma_x_slice_V/1e6);
           xlabel('y'); ylabel('\sigma_x [MPa]');xlim([-W,0])
Integral_Moment_x =  W*trapz(dist,dist.*sigma_x_slice_V)
Integral_tau_xy   = W*trapz(dist,tau_xy_slice_V)

vonMises = EvaluateVonMises(sigma_x,sigma_y,tau_xy);

figure(5); FEMtrisurf(FEMmesh,vonMises/1e6);
           xlabel('x'); ylabel('y'); zlabel('von Mises [MPa]'); view(160,25)
            colorbar(); shading interp

figure(6); clf
           FEMtricontour(FEMmesh,vonMises/1e6,1e1*[0:0.5:6]);
	   xlabel('x'); ylabel('y'); title('von Mises [MPa]');
           caxis(1e2*[0 0.7]); axis equal; colorbar()
           hold on; plot([Domain(:,1);Domain(1,1)],[Domain(:,2);Domain(1,2)],'color','black','linewidth',1); hold off

dist = linspace(-W,gap,100)'; HH = H/2*ones(size(dist));
vonMises_slice = FEMgriddata(FEMmesh,vonMises,dist,dist);
figure(7); plot(dist,vonMises_slice*1e-6); xlabel('x,y'); ylabel('von Mises [MPa]')

