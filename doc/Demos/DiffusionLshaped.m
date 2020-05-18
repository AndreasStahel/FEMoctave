nodes = [0,0,-2;1,0,-2;1,1,-2;-1,1,-2;-1,-1,-2;0,-1,-2];
FEMmesh = CreateMeshTriangle('Ldomain',nodes,0.02);
%%FEMmesh = MeshUpgrade(FEMmesh);  %% uncomment to use quadratic elements

u = BVP2Dsym(FEMmesh,1,0,1,0,0,-2);

figure(1);
FEMtrimesh(FEMmesh.elem,FEMmesh.nodes(:,1),FEMmesh.nodes(:,2),u);
xlabel('x'); ylabel('y'); view(-30,30); grid on
