%% script file to test code
clear *
N = 10;  Triangle = 1

%% Dirichlet and Neumann BC
if Triangle
  FEMmesh = CreateMeshTriangle('test1',[0 0 -2;1 0 -1; 1 1 -1; 0 1 -2],0.75/N^2);
  FEMmesh = MeshUpgrade(FEMmesh);
%%  FEMmesh1 = CreateMeshTriangle('test1',[0 0 -2;1 0 -1; 1 1 -1; 0 1 -2],0.75/N^2/4);
  FEMmesh1 = MeshQuad2Linear(FEMmesh);
  nDOFTri = [FEMmesh.nDOF, FEMmesh1.nDOF]
else
  FEMmesh = CreateMeshRect(linspace(0,1,N+1),linspace(0,1,N+1),-2,-1,-2,-1);
  FEMmesh = MeshUpgrade(FEMmesh);
%%  FEMmesh1 = CreateMeshRect(linspace(0,1,2*N+1),linspace(0,1,2*N+1),-2,-1,-2,-1);
  FEMmesh1 = MeshQuad2Linear(FEMmesh);
  nDOFRect = [FEMmesh.nDOF, FEMmesh1.nDOF]
endif


function y = uSol(xy)
  y = besselj(0,sqrt(xy(:,1).^2+xy(:,2).^2));
endfunction
%% the Besselj(0,r) function is an eigenfunction with eigenvalue 1
function y = LapuSol(xy)
  y = 2*besselj(0,sqrt(xy(:,1).^2+xy(:,2).^2));
endfunction

function y = aa(xy)
  y = ones(length(xy),1);
endfunction

u     = BVP2Dsym(FEMmesh ,"aa",1,"LapuSol","uSol",0,0);
u_lin = BVP2Dsym(FEMmesh1,"aa",1,"LapuSol","uSol",0,0);

uExact = uSol(FEMmesh.nodes); uExact1 = uSol(FEMmesh1.nodes);

figure(1);
FEMtrimesh(FEMmesh,uExact);
xlabel("x"); ylabel("y"); title("exact solution"); view([45,25])


figure(2);
FEMtrimesh(FEMmesh,u-uExact);
title('error, quadratic elements');xlabel("x");ylabel("y")
view([45,25])

[xi,yi] = meshgrid(linspace(0,1,100));
ui_exact = besselj(0,sqrt(xi.^2+yi.^2));
ui = FEMgriddata(FEMmesh,u,xi,yi);

figure(3)
mesh(xi,yi,ui-ui_exact)
xlabel('y'); ylabel('y'); title('difference')

figure(12);
FEMtrimesh(FEMmesh1,u_lin-uExact1);
title('error, linear elements');xlabel("x");ylabel("y")

MaxMeanL2Error2A = [max(abs(u-uExact)),mean(abs(u-uExact)),sqrt(mean((u-uExact).^2))]
MaxMeanL2Error1A = [max(abs(u_lin-uExact1)),mean(abs(u_lin-uExact1)),sqrt(mean((u_lin-uExact1).^2))]

x = FEMmesh.nodes(:,1); y = FEMmesh.nodes(:,2);
r = sqrt(x.^2+y.^2); phi = mod(atan2(y,x),2*pi);
ux_exact = -cos(phi).*besselj(1,r);
uy_exact = -sin(phi).*besselj(1,r);

x1 = FEMmesh1.nodes(:,1); y1 = FEMmesh1.nodes(:,2);
r1 = sqrt(x1.^2+y1.^2); phi1 = mod(atan2(y1,x1),2*pi);
ux1_exact = -cos(phi1).*besselj(1,r1);
uy1_exact = -sin(phi1).*besselj(1,r1);

[ux,uy]   = FEMEvaluateGradient(FEMmesh,uExact);
[ux1,uy1] = FEMEvaluateGradient(FEMmesh1,uExact1);

figure(13);
FEMtrimesh(FEMmesh,ux_exact);
xlabel("x"); ylabel("y"); title('u_x exact')
view([30,30])


figure(23);
FEMtrimesh(FEMmesh,ux-ux_exact);
xlabel("x"); ylabel("y"); title('error u_x, quadratic elements')
view([45,25])

figure(33);
FEMtrimesh(FEMmesh1.ux1-ux1_exact);
xlabel("x"); ylabel("y"); title('error u_x, linear elements')
view([30,30])
view([45,25])

