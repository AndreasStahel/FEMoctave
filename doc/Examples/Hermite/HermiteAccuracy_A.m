N = 11;
Mesh = CreateMeshRect(linspace(0,1,N),linspace(0,1,N),-1,-2,-1,-2);
MeshHermite = MeshUpgrade(Mesh,'Hermite');
MeshCubic = MeshUpgrade(Mesh,'cubic'); MeshCubic = MeshCubic2Linear(MeshCubic);
MeshCubic = MeshUpgrade(MeshCubic,'cubic');  %% a rather fine mesh with cubic elements

a = 1; b0 = 0; bx = 0; by = 0; gD = 1; gN1 = 0; gN2 = 0;
function res = f(xy,dummy)
   res = 1*xy(:,1).^2;
endfunction

u_ref  = BVP2D(MeshCubic,  a,b0,bx,by,"f",gD,gN1,gN2);
u_Herm = BVP2D(MeshHermite,a,b0,bx,by,"f",gD,gN1,gN2);
figure(1); FEMtrimesh(MeshHermite,u_Herm); xlabel('x'); ylabel('y'); zlabel('u_{Herm}')
%% use griddata on small part of the domain
[xx,yy] = meshgrid(linspace(0,0.2,2*N),linspace(0,0.1,2*N));
u_ref_g  = FEMgriddata(MeshCubic,u_ref,xx,yy);
u_Herm_g = FEMgriddata(MeshHermite,u_Herm,xx,yy);
figure(2); mesh(xx,yy,u_Herm_g); xlabel('x'); ylabel('y'); zlabel('u_{Herm}')
Difference_RMS = sqrt(mean((u_Herm_g(:)-u_ref_g(:)).^2))
figure(21); mesh(xx,yy,u_Herm_g-u_ref_g); xlabel('x'); ylabel('y'); zlabel('u_{Herm}-u_{ref}')
levels = 1+linspace(-0.001,0.01,23);
figure(22); contour(xx,yy,u_ref_g,levels);  xlabel('x'); ylabel('y'); title('u_{ref}')
figure(23); contour(xx,yy,u_Herm_g,levels); xlabel('x'); ylabel('y'); title('u_{Herm}')
%%  along the two edges
N1 = 1001;
x = linspace(0,1,N1); y = zeros(size(x));
[u_Herm_gx,ux,uy] = FEMgriddata(MeshHermite,u_Herm,x,y);
figure(31); plot(x,u_Herm_gx); xlabel('x'); ylabel('u_{Herm}(x,0)')
figure(32); plot(x,ux); xlabel('x'); ylabel('du/dx(x,0)')
[u_Herm_gy,ux,uy] = FEMgriddata(MeshHermite,u_Herm,y,x);
figure(33); plot(x,u_Herm_gy); xlabel('y'); ylabel('u_{Herm}(0,y)')
figure(34); plot(x,uy); xlabel('y'); ylabel('du/dy(0,y)')
