%% define the domain
clear *
R = 2; R_in = 1.0;
area = 0.001;
z = linspace(-1+sqrt(area),1-sqrt(area),21)'; r = R-R_in*cos(pi/2*z).^2; b = -2*ones(size(z));
zr = [-2 0 -1; -2 R -2; -1 R -2; [z,r,b]; 1 R -2; 2 R -1; 2 0 -2];
if 0      %% linear elements
  FEMmesh = CreateMeshTriangle('PotentialFlow',zr,area);
elseif 1 %% quadratic elements
  FEMmesh = CreateMeshTriangle('PotentialFlow',zr,4*area);
  FEMmesh = MeshUpgrade(FEMmesh,'quadratic');
else    %% cubic elements
  FEMmesh = CreateMeshTriangle('PotentialFlow',zr,9*area);
  FEMmesh = MeshUpgrade(FEMmesh,'cubic');
endif

z = FEMmesh.nodes(:,1);  z = FEMmesh.nodes(:,2);
function res = gD(zr)       res = -zr(:,1)/2;  endfunction
function res = a_coeff(zr)  res = zr(:,2);     endfunction


u = BVP2Dsym(FEMmesh,'a_coeff',0,0,'gD',0,0);
figure(1); FEMtrimesh(FEMmesh,u)
           xlabel('z'); ylabel('r'); zlabel('potential')
           view(20,40)

[zz,rr] = meshgrid(linspace(-2,2-0.01,35),linspace(0,R-0.01,41));
[u_int,uz_int,ur_int] = FEMgriddata(FEMmesh,-u, zz, rr);

figure(2); quiver(zz,rr,uz_int,ur_int)
           xlabel('z'); ylabel('r');
           hold on; plot([zr(:,1);-2],[zr(:,2);0],'k'); hold off
           xlim([-2,2]); ylim([0,R]);

[uz,ur] = FEMEvaluateGradient(FEMmesh,u);
figure(4); FEMtrimesh(FEMmesh,sqrt(uz.^2+ ur.^2))
           xlabel('z'); ylabel('r'); zlabel('v=|grad u|')
           zlim([0 1]); caxis([0,1])


figure(14); mesh(zz,rr,sqrt(uz_int.^2+ ur_int.^2))
            xlabel('z'); ylabel('r'); zlabel('v=|grad u|')
            zlim([0 1]); caxis([0 1])

zz = linspace(-2,2,101); rr = 0.5*ones(101,1);
[u_int,uz_int,ur_int] = FEMgriddata(FEMmesh,-u,zz,rr);

figure(3); plot(zz,uz_int)
           xlabel('z'); ylabel('horizontal velocity');
           ylim([0 1.1*max(uz_int)])

figure(5)
FEMtricontour(FEMmesh,sqrt(uz.^2+ ur.^2),31)
xlabel('z'); ylabel('r'); zlabel('|grad u|')
hold on; plot([zr(:,1);-2],[zr(:,2);0],'k'); hold off
xlim([-2 2]); ylim([0 R]);
%%colorbar
axis equal

rr = linspace(0,R); zz = -1.9*ones(size(rr));
vz = FEMgriddata(FEMmesh,-uz, zz, rr);
Flux_inlet = trapz(rr,rr.*vz)*2*pi
rr = linspace(0,R-R_in); zz = 0*ones(size(rr));
vz = FEMgriddata(FEMmesh,-uz,zz,rr);
Flux_middle = trapz(rr,rr.*vz)*2*pi
rr = linspace(0,R); zz = 1.9*ones(size(rr));
vz = FEMgriddata(FEMmesh,-uz, zz, rr);
Flux_outlet = trapz(rr,rr.*vz)*2*pi

