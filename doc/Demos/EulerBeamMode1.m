clear *
L = 0.20; H = 0.01; W = 0.01; rho = 2.7e3;
E = 70e9; nu = 0.33; %% Aluminum
I2 = 1/12*H^3*W;

f = @(z) 1+cos(z).*cosh(z);  %% clamped at x=0, free at x=L
z0 = fsolve(f,pi/2);
freqEuler = z0^2*sqrt(E*I2/(rho*H*W))/(2*pi*L^2)
Nx = 20; Ny = 2;
Mesh = CreateMeshRect(linspace(0,L,Nx+1),linspace(0,+H,Ny+1),-22,-22,-11,-22);
Mesh = MeshUpgrade(Mesh,'quadratic');
[la,u1,u2] = PlaneStressEig(Mesh,E,nu,rho,1);
freqFEM = sqrt(la)/(2*pi)
u1 = u1/max(abs(u2))/100; u2 = u2/max(abs(u2))/100;
figure(1);FEMtrimesh(Mesh,u1); xlabel('x'); ylabel('y'); zlabel('u_1')
figure(2);FEMtrimesh(Mesh,u2); xlabel('x'); ylabel('y'); zlabel('u_2')
