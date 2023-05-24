%% solve a static heat problem in a ball, using the radius r only as variable
R = 3; W = 0.1; N = 10;
FEMmesh = CreateMeshRect(linspace(1,R,N+1),[-W/2,+W/2],-2,-2,-1,-1);
%%FEMmesh = MeshUpgrade(FEMmesh,'quadratic');
%%FEMmesh = MeshUpgrade(FEMmesh,'cubic');

function res = a(r) ;  res = r(:,1).^2;      endfunction
function res = gD(r);  res = sign(r(:,1)-1); endfunction

u3D = BVP2Dsym(FEMmesh,'a',0,0,'gD',0,0);
figure(1); FEMtrimesh(FEMmesh,u3D)
           xlabel('r'); ylabel('dummy'); zlabel('temperature u')

r = linspace(1,R);
u = FEMgriddata(FEMmesh,u3D,r,0*r);
u_exact = R/(R-1)*(r-1)./r;

figure(2); plot(r,u,r,u_exact)
            xlabel('radius r'); ylabel('temperature u')
            legend('u_{FEM}','u_{exact}','location','northwest')

figure(3); plot(r,u-u_exact)
            xlabel('radius r'); ylabel('u_{FEM}-u_{exact}')
