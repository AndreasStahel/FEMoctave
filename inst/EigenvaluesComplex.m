## -*- texinfo -*-
## @deftypefn  {} {} EigenvaluesComplex.m ()
##
## This is a demo file  inside the `doc/Demos/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

FEMmesh = CreateMeshRect(linspace(0,pi,21),linspace(0,pi,21),-1,-1,-1,-1);
FEMmesh = MeshUpgrade(FEMmesh,'quadratic');

%%%%%%% solve the eigenvalue problem, show the eigenvalues
[la,ve] = BVP2Deig(FEMmesh,1,i,1,0,7,"type","complex");
eigenvalues = la
figure(1); FEMtrimesh(FEMmesh,real(ve(:,3)));  xlabel("x"); ylabel("y"); zlabel("real(u)")
figure(2); FEMtrimesh(FEMmesh,imag(ve(:,3)));  xlabel("x"); ylabel("y"); zlabel("imag(u)")

