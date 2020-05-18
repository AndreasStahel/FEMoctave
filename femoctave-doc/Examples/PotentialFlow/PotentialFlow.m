%% define the domain
xy =[0 0 -2;5 0 -1;5 2 -2;3 2 -2;3 0.5 -2;2 0.5 -2;1 2 -2;0 2 -1];
if 0  %% linear elements
  FEMmesh = CreateMeshTriangle('PotentialFlow',xy,0.003);
else  %% quadratic elements
  FEMmesh = CreateMeshTriangle('PotentialFlow',xy,4*0.003);
  FEMmesh = MeshUpgrade(FEMmesh);
endif

x = FEMmesh.nodes(:,1);  y = FEMmesh.nodes(:,2);
function res = gD(xy)   res = 1-xy(:,1)/5; endfunction
u = BVP2Dsym(FEMmesh,1,0,0,'gD',0,0);
figure(1)
trimesh(FEMmesh.elem,FEMmesh.nodes(:,1),FEMmesh.nodes(:,2),u)
xlabel('x'); ylabel('y'); zlabel('potential')

[xx,yy] = meshgrid(linspace(0,5-0.01,25),linspace(0,2-0.01,21));
[u_int,ux_int,uy_int] = FEMgriddata(FEMmesh,-u, xx, yy);

figure(2)
quiver(xx,yy,ux_int,uy_int)
xlabel('x'); ylabel('y');
hold on; plot([xy(:,1);0],[xy(:,2);0],'k'); hold off; axis equal


[ux,uy] = FEMEvaluateGradient(FEMmesh,u);
figure(4)
FEMtrimesh(FEMmesh.elem,FEMmesh.nodes(:,1),FEMmesh.nodes(:,2),sqrt(ux.^2+ uy.^2))
xlabel('x'); ylabel('y'); zlabel('v=|grad u|'); zlim([0 0.5])

figure(14)
mesh(xx,yy,sqrt(ux_int.^2+ uy_int.^2))
xlabel('x'); ylabel('y'); zlabel('v=|grad u|'); zlim([0 0.5])

xx = linspace(0,5,101); yy = 0.25*ones(101,1);
[u_int,ux_int,uy_int] = FEMgriddata(FEMmesh,-u,xx,yy);
figure(3)
plot(xx,ux_int)
xlabel('x'); ylabel('horizontal velocity'); ylim([0 0.5])

figure(5)
if size(FEMmesh.elem,2)==6 FEMmesh = MeshQuad2Linear(FEMmesh); endif
tricontour(FEMmesh.elem,FEMmesh.nodes(:,1),FEMmesh.nodes(:,2),sqrt(ux.^2+ uy.^2),21)
xlabel('x'); ylabel('y'); zlabel('| grad u|')
hold on; plot([xy(:,1);0],[xy(:,2);0],'k'); hold off
xlim([0 5]); ylim([0 2]); axis equal

yy = linspace(0,2); xx = zeros(size(yy));
vx = FEMgriddata(FEMmesh,-ux, xx, yy);
Flux_inlet_ = trapz(yy,vx)
yy = linspace(0,0.5); xx = 2.5*ones(size(yy));
vx = FEMgriddata(FEMmesh,-ux,xx,yy);
Flux_middle = trapz(yy,vx)
yy = linspace(0,2); xx = 5*ones(size(yy));
vx = FEMgriddata(FEMmesh,-ux, xx, yy);
Flux_outlet = trapz(yy,vx)

