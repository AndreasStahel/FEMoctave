pkg load femoctave
R0 = 1; R1 = 3.5; N = 20;
%% create the mesh
FEMmesh = CreateMeshRect(linspace(R0,R1,N),linspace(0,pi,N),-2,-2,-2,-1);
figure(1); FEMtrimesh(FEMmesh); xlabel('radius r'); ylabel('angle \phi')
function xy = Deform(r_phi)
  r = r_phi(:,1); phi = r_phi(:,2);
  xy = [r.*cos(phi),r.*sin(phi)];
endfunction
figure(11); FEMtrimesh(MeshDeform(FEMmesh,'Deform')); axis equal; xlabel('x'); ylabel('y');
FEMmesh = MeshUpgrade(FEMmesh,"cubic");

%% define the functions
function res = gD(r_phi)
  res = -r_phi(:,1).*cos(r_phi(:,2));
endfunction

function res = a(r_phi)
   r = r_phi(:,1);
   res = [r,1./r,0*r];
endfunction

%% solve the BVP
Phi = BVP2Dsym(FEMmesh,"a",0,0,"gD",0,0);

figure(2); FEMtrimesh(FEMmesh,Phi)
           xlabel('r'); ylabel('\phi'); zlabel('potential \Phi')
[u_r,u_phi] = FEMEvaluateGradient(FEMmesh,Phi);
v = sqrt(u_r.^2+(u_phi./FEMmesh.nodes(:,1)).^2);
figure(3); FEMtrimesh(FEMmesh,v)
           xlabel('r'); ylabel('\phi'); zlabel('velocity');
FEMmeshCartesian = MeshDeform(MeshCubic2Linear(FEMmesh),'Deform');
figure(33); FEMtrimesh(FEMmeshCartesian,v)
            xlabel('x'); ylabel('y'); zlabel('speed');
            xlim([-R1,+R1]);ylim([0,R1]); zlim([0 2]);
figure(34); FEMtrimesh(FEMmeshCartesian,v); view([0,90])
            xlabel('x'); ylabel('y'); colorbar; axis equal; axis equal

%% evaluate the  gradient as vector field
r = FEMmesh.nodes(:,1); phi = FEMmesh.nodes(:,2);
x = r.*cos(phi); y = r.*sin(phi);
u_x = +cos(phi).*u_r - sin(phi).*u_phi./r;
u_y = +sin(phi).*u_r + cos(phi).*u_phi./r;
figure(4); quiver(x,y,-u_x,-u_y)
           xlabel('x'); ylabel('y'); xlim([-2,2]); ylim([0,2]); axis equal

%% find the flow lines
[xx,yy] = meshgrid(linspace(-2,2,51),linspace(0,3,25));
r_int = sqrt(xx.^2+yy.^2); phi_int = atan2(yy,xx);
[u_int,ur_int,uphi_int] = FEMgriddata(FEMmesh,-Phi,r_int,phi_int);
ux_int = +cos(phi_int).*ur_int - sin(phi_int).*uphi_int./r_int;
uy_int = +sin(phi_int).*ur_int + cos(phi_int).*uphi_int./r_int;
ux_int(find(isnan(ux_int))) = 0;  uy_int(find(isnan(uy_int))) = 0;

figure(5); clf; NN = 50; alpha = linspace(0,pi);
         streamline(xx,yy,ux_int,uy_int,-2*ones(1,NN),...
                    linspace(0.001,3,NN),[0.1,10000]);
         hold on; plot(cos(alpha),sin(alpha),'k'); hold off; axis equal
         xlabel('x'); ylabel('y'); xlim([-2 2]); ylim([0,2.5]);
