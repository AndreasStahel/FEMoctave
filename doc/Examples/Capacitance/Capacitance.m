FEMmesh = ReadMeshTriangle('capacitance.1');
%%FEMmesh = MeshUpgrade(FEMmesh);
x = FEMmesh.nodes(:,1); y = FEMmesh.nodes(:,2); 

figure(1)
FEMtrimesh(FEMmesh)

function res = a(xy)     res = xy(:,1);      endfunction
function res = Volt(xy)  res = xy(:,2)>0.1;  endfunction

u = BVP2Dsym(FEMmesh,'a',0,0,'Volt',0,0);
figure(2)
FEMtrimesh(FEMmesh,u);
view([38,48])
xlabel('radius r'); ylabel('height z'); zlabel('voltage')

if size(FEMmesh.elem,2)==6
  MeshTri = MeshQuad2Linear(FEMmesh);
else
  MeshTri = FEMmesh;
endif
figure(3)
tricontour(MeshTri.elem,MeshTri.nodes(:,1),MeshTri.nodes(:,2),u,21);
xlabel('radius r'); ylabel('height z');

[ux,uy] = FEMEvaluateGradient(FEMmesh,u);
xi = linspace(0,2.5,101)'; yi = 0.0*ones(101,1);
uy_i = FEMgriddata(FEMmesh,uy,xi,yi);

figure(4)
plot(xi,uy_i)
xlabel('radius r'); ylabel('u_z'); ylim([-1,6])

Integral = [2*pi*trapz(xi,xi.*uy_i), pi*1^2/0.2]
