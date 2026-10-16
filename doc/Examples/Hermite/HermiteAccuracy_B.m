## -*- texinfo -*-
## @deftypefn  {} {} HermiteAccuracy_B.m
##
## This is a demo file inside the `doc/Examples/Hermite/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

%% solution is a bell surface on a full circle
N = 20; R = 2;
angles = linspace(0,2*pi,2*N)'; angles = angles(1:end-1);
xy = [R*cos(angles),R*sin(angles),-ones(size(angles))];
Mesh = CreateMeshTriangle('circ',xy,4*R^2/N^2);
x = Mesh.nodes(:,1); y = Mesh.nodes(:,2);
a = 1; b0 = 0; bx = 0; by = 0; gN1 = 0; gN2 = 0;
function res = gD(xy)
  res = exp(-sum(xy.^2,2));
endfunction

function res = f(xy,d)
  x = xy(:,1); y = xy(:,2);
  res = -(4*x.^2+4*y.^2-4).*exp(-sum(xy.^2,2));
endfunction

tic();
Mesh3 = MeshUpgrade(Mesh,'cubic');
u3 = BVP2D(Mesh3,a,b0,bx,by,'f','gD', gN1,gN2);
Time_cubic = toc()
figure(1); FEMtrimesh(Mesh3,u3); xlabel('x'); ylabel('y'); zlabel('u_3')
x3 = Mesh3.nodes(:,1); y3 = Mesh3.nodes(:,2);
Difference_u3 = mean(abs(u3-gD([x3,y3])))
%% Hermite elements
tic();
MeshH = MeshUpgrade(Mesh,'Hermite');
uH = BVP2D(MeshH,a,b0,bx,by,'f','gD', gN1,gN2);
Time_Hermite = toc()
figure(2); FEMtrimesh(MeshH,uH);
           xlabel('x'); ylabel('y'); zlabel('u_{H}')
tic();
uHI = FEMgriddata(MeshH,uH,x,y);
Time_HermInterpol = toc()
figure(3); FEMtrimesh(Mesh,uHI);
           xlabel('x'); ylabel('y'); zlabel('u_{HI}')
Difference_uHI = mean(abs(uHI-gD([x,y])))
%% interpolation for Hermite
tic();
uH3 = FEMgriddata(MeshH,uH,x3,y3);
Time_HermInterpol3 = toc()
figure(4); FEMtrimesh(Mesh3,uH3);
           xlabel('x'); ylabel('y'); zlabel('u_{H3}')
Difference_uH3 = mean(abs(uH3-gD([x3,y3])))
figure(5); FEMtrimesh(Mesh3,uH3-u3);
           xlabel('x'); ylabel('y'); zlabel('difference')
