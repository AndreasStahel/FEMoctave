R = 0.1;
if 0  %% nonuniform mesh
  Mesh = CreateMeshTriangle('AxiSymm',[0 0 -21; R 0 -32; R 2*R -22; 0 2*R -12],1e-4);
else
  Mesh = CreateMeshRect(linspace(0,R,10),linspace(0,2*R,20),-21,-22,-12,-32);
endif
Mesh = MeshUpgrade(Mesh,'quadratic');

function res = force(rz)
  R = 0.1; P = 1e6;
  res = -P*max(R^2-rz(:,2).^2,0);
endfunction


E = 110e9; nu = 0.35; f = {0,0}; gD = {0,0}; gN = {'force',0};
[ur,uz] = AxiStress(Mesh,E,nu,f,gD,gN);

factor = 0.1*R/max(sqrt(ur.^2+uz.^2));
figure(1); ShowDeformation(Mesh,ur,uz,factor); xlabel('r'); ylabel('z'); axis equal

[sigma_x,sigma_y,sigma_z,tau_xz] = EvaluateStressAxi(Mesh,ur,uz,E,nu);
figure(12); FEMtrimesh(Mesh,sigma_x)
            xlabel('r'); ylabel('z'); zlabel('\sigma_x')
figure(13); FEMtrimesh(Mesh,sigma_y)
            xlabel('r'); ylabel('z'); zlabel('\sigma_y')
figure(14); FEMtrimesh(Mesh,sigma_z)
            xlabel('r'); ylabel('z'); zlabel('\sigma_z')


vonMises = EvaluateVonMises(sigma_x,sigma_y,sigma_z,tau_xz);
figure(15); FEMtrimesh(Mesh,vonMises)
            xlabel('r'); ylabel('z'); zlabel('von Mises'); view([-125,30])

