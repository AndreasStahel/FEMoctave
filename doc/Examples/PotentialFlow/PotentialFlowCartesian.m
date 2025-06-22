pkg load femoctave
R0 = 1; R1 = 3.5; N = 20;
%% create the mesh
FEMmesh = CreateMeshRect(linspace(R0,R1,N),linspace(0,pi,N),-2,-2,-2,-1);
function xy = Deform(r_phi)
  r = r_phi(:,1); phi = r_phi(:,2);
  xy = [r.*cos(phi),r.*sin(phi)];
endfunction
FEMmesh = MeshDeform(FEMmesh,'Deform');
figure(1); FEMtrimesh(FEMmesh); axis equal; xlabel('x'); ylabel('y');
FEMmesh = MeshUpgrade(FEMmesh,"cubic");
%%FEMmesh = MeshUpgrade(FEMmesh,"quadratic");

%% define the function
function res = gD(xy) res = -xy(:,1); endfunction

%% solve the BVP
Phi = BVP2Dsym(FEMmesh,1,0,0,"gD",0,0);

figure(2); FEMtrimesh(FEMmesh,Phi)
           xlabel('r'); ylabel('\phi'); zlabel('potential \Phi')
[u_x,u_y] = FEMEvaluateGradient(FEMmesh,Phi);
v = sqrt(u_x.^2+u_y.^2);
figure(3); FEMtrimesh(FEMmesh,v)
           xlabel('x'); ylabel('y'); zlabel('velocity');
figure(34); FEMtrimesh(FEMmesh,v); view([0,90])
            xlabel('x'); ylabel('y'); colorbar; axis equal

%% evaluate the  gradient as vector field
x = FEMmesh.nodes(:,1); y = FEMmesh.nodes(:,2);
figure(4); quiver(x,y,-u_x,-u_y)
           xlabel('x'); ylabel('y'); xlim([-2,2]); ylim([0,2]); axis equal

%% find the flow lines
[xx,yy] = meshgrid(linspace(-2,2,51),linspace(0,3,25));
[u_int,ux_int,uy_int] = FEMgriddata(FEMmesh,-Phi,xx,yy);
ux_int(find(isnan(ux_int))) = 0;  uy_int(find(isnan(uy_int))) = 0;

figure(5); clf; NN = 50; alpha = linspace(0,pi);
         streamline(xx,yy,ux_int,uy_int,-2*ones(1,NN),...
                    linspace(0.001,3,NN),[0.1,10000]);
         hold on; plot(cos(alpha),sin(alpha),'k'); hold off; axis equal
         xlabel('x'); ylabel('y'); xlim([-2 2]); ylim([0,2.5]);
