
Mesh = CreateMeshTriangle('square',[-1,-1,-2;1,-1,-2;1 1 -2; -1 1 -2],0.05);
figure(1); FEMtrimesh(Mesh);
Mesh = MeshUpgrade(Mesh,'quadratic');
function res = f(xy)
  res = sin(pi*xy(:,1));
endfunction
u1 = BVP2D(Mesh,1,0,0,0,'f',7,0,0);
MeanSolution1 = mean(u1)

Mesh = MeshAddConstraint(Mesh,[0,0],-1);
u2 = BVP2D(Mesh,1,0,0,0,'f',7,0,0);
figure(1); FEMtrimesh(Mesh,u2); xlabel('x'); ylabel('y'); zlabel('u')
MeanSolution2 = mean(u2)
MeanSolution2Integration = FEMIntegrate(Mesh,u2)/4
