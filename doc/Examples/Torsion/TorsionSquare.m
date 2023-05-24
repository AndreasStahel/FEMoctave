N = 10;
l = sqrt(pi)/2;  al = 1; %%al = sqrt(2);
if 0 %% uniform mesh
  Mesh = CreateMeshRect(al*linspace(-l,l,N),1/al*linspace(-l,l,N),-1,-1,-1,-1);
else %% non uniform mesh
  Mesh = CreateMeshTriangle('Torsion',
  [-al*l -1/al*l -1; al*l -1/al*l -1; al*l 1/al*l -1; -al*l 1/al*l -1],pi/2/N^2);
endif
Mesh = MeshUpgrade(Mesh);

chi = BVP2Dsym(Mesh,1,0,2,0,0,0);
[chiGP,gradChi] = FEMEvaluateGP(Mesh,chi);
xGP = Mesh.GP(:,1); yGP = Mesh.GP(:,2);
f = xGP.*gradChi(:,1) + yGP.*gradChi(:,2);
J = FEMIntegrate(Mesh,-f)

[chi_x,chi_y] = FEMEvaluateGradient(Mesh,chi);
Stress = sqrt(chi_x.^2 + chi_y.^2);
figure(1)
FEMtrisurf(Mesh,Stress)
xlabel('x'); ylabel('y');

MaxStress = max(Stress)
%print -dpdfwrite doc/TorsionSquare.pdf
%print -dpdfwrite doc/TorsionRectangle.pdf

