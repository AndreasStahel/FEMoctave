FEMmesh = CreateMeshRect([0:0.1:1],[0:0.1:2],-1,-2,-1,-2);
FEMmesh = MeshUpgrade(FEMmesh);  %% uncomment to use quadratic elements

u = BVP2Dsym(FEMmesh,1,0,0.25,0,0,0);

figure(1)
FEMtrimesh(FEMmesh.elem,FEMmesh.nodes(:,1),FEMmesh.nodes(:,2),u);
xlabel('x'); ylabel('y');    
