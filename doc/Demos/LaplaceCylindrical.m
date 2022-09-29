FEMmesh = CreateMeshRect(linspace(0,2,20),linspace(-1,2,30),-1,-1,-2,-2);
%%FEMMesh = CreateMeshTriangle('test',[0 -1 -1; 2 -1 -2; 2 2 -1; 0 2 -2],0.1); 
FEMmesh = MeshUpgrade(FEMmesh, 'quadratic');  %% uncomment to use quadratic elements

function res = f(rz)   res = rz(:,1)*2.*rz(:,2); endfunction
function res = b0(rz)  res = 10*rz(:,1);         endfunction
function res = a(rz)   res = rz(:,1);            endfunction
function res = g2(rz)  res = -1*rz(:,1)/2;       endfunction

u = BVP2Dsym(FEMmesh,'a','b0','f',0,'g2',0);

FEMtrimesh(FEMmesh,u);
xlabel('\rho'); ylabel('z');
