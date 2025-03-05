%% Poisseuille flow in an annulus, 2D solution
G = 1; mu = 1; R1 = 1; R2 = 2;
angles = linspace(0,2*pi,201)'; angles = angles(1:end-1);
xy = [R2*cos(angles),R2*sin(angles),-ones(size(angles))];
Hole.name   = 'Hole';
Hole.border = [R1*cos(angles),R1*sin(angles),-ones(size(angles))];
Hole.point  = [0,0];
Mesh = CreateMeshTriangle('Annulus',xy,0.01,Hole);
Mesh = MeshUpgrade(Mesh,'quadratic');

u = BVP2D(Mesh,1,0,0,0,G/mu,0,0,0);

figure(1); FEMtrimesh(Mesh,u)
           xlabel('x'); ylabel('y'); zlabel('velocity u')
Flow = FEMIntegrate(Mesh,u)
