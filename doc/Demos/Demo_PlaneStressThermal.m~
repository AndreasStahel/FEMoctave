E = 110e3; nu = 0.35; alpha = 1.7e-5; %% copper, length in mm
L = 0.1;  N = 11;

Mesh = CreateMeshRect(linspace(0,L,N),linspace(0,L,N),-11,-11,-11,-11);
figure(1); FEMtrimesh(Mesh)
% Mesh = MeshUpgrade(Mesh,'quadratic');
Mesh = MeshUpgrade(Mesh,'cubic');

g1 = @(xy)(0.005*xy(:,1).^2); g2 = @(xy)(0.005*xy(:,2).^2);
alphaDeltaT = alpha*30;  %% approx 0.0005

[u1,u2] = PlaneStress(Mesh,E,nu,{0,0},{g1,g2},{0,0},'thermal',alphaDeltaT);
figure(2); FEMtrimesh(Mesh,u1); xlabel('x'); ylabel('y'); zlabel('u_1')
figure(3); FEMtrimesh(Mesh,u2); xlabel('x'); ylabel('y'); zlabel('u_2')

[eps_xx,eps_yy,eps_xy] = EvaluateStrain(Mesh,u1,u2);
figure(4); FEMtrimesh(Mesh,eps_xx); xlabel('x'); ylabel('y'); zlabel('\epsilon_{xx}')
figure(5); FEMtrimesh(Mesh,eps_yy); xlabel('x'); ylabel('y'); zlabel('\epsilon_{yy}')

[sigma_x,sigma_y,tau_xy] = EvaluateStress(Mesh,u1,u2,E,nu,'thermal',alphaDeltaT);
figure(6); FEMtrimesh(Mesh,sigma_x); xlabel('x'); ylabel('y'); zlabel('\sigma_x')
figure(7); FEMtrimesh(Mesh,sigma_y); xlabel('x'); ylabel('y'); zlabel('\sigma_y')
figure(8); FEMtricontour(Mesh,sigma_x,41); xlabel('x'); ylabel('y');
           colorbar(); axis equal
figure(9); FEMtricontour(Mesh,sigma_y,41); xlabel('x'); ylabel('y');
           colorbar(); axis equal

MaxSigma_x  = max(abs(sigma_x))
MeanSigma_x = mean(sigma_x)

