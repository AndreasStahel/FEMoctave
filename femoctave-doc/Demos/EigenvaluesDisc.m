clear *
%% create a circle with mesh
xM = 0; yM = 0; R = 1; N = 1*40;
alpha = linspace(0,N/(N+1)*2*pi,N)';
xy = [xM+R*cos(alpha),yM+R*sin(alpha),-ones(size(alpha))];

FEMmesh = CreateMeshTriangle("circle",xy,0.001);
%%FEMmesh = MeshUpgrade(FEMmesh);

%%%%%%% solve the eigenvalue problem, show the eigenvalues
if 1
  [la,ve,errorbound] = BVP2Deig(FEMmesh,1,0,1,0,4,1e-3);
  eigenvalues = la
  errorbound
  figure(1);
  FEMtrimesh(FEMmesh.elem,FEMmesh.nodes(:,1),FEMmesh.nodes(:,2),ve(:,4));
  xlabel("x"); ylabel("y");
else
  la = BVP2Deig(FEMmesh,1,0,1,0,4,1e-8)
endif

exact_values = [fsolve(@(x)besselj(0,x),2.3),fsolve(@(x)besselj(1,x),3.8),fsolve(@(x)besselj(2,x),5)].^2
