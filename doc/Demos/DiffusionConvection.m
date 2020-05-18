clear *
if 0  %% linear elements
  FEMmesh = CreateMeshRect(linspace(0,2,51),linspace(0,2,51),-1,-1,-1,-1);
  u = BVP2D(FEMmesh,1,0,10,5,0.1,0,0,0);
  FEMmesh.nDOF
else  %% quadratic elements
  FEMmesh = CreateMeshRect(linspace(0,2,26),linspace(0,2,26),-1,-1,-1,-1);
  FEMmesh = MeshUpgrade(FEMmesh);   %% make a mesh with elements of order 2
  FEMmesh.nDOF
  u = BVP2D(FEMmesh,1,0,10,5,0.1,0,0,0);
  FEMmesh = MeshQuad2Linear(FEMmesh);  %% convert to identical mesh with elements of order 1
  FEMmesh.nDOF
endif

figure(1)
tricontour(FEMmesh.elem,FEMmesh.nodes(:,1),FEMmesh.nodes(:,2),u,10);
colorbar()
xlabel('x'); ylabel('y'); grid on
