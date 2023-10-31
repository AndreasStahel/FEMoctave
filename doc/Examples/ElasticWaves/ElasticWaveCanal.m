x = linspace(0,1,21)'; y = cos(pi/2*x).^2; bc = -22*ones(size(x));
MeshData = [-1,0.5,-22;x,0.5*y,bc;2.5,0,-11;2.5,-1,-22;-1,-1,-22];

x = MeshData(:,1); y = MeshData(:,2);
Mesh = CreateMeshTriangle('Canal',MeshData,0.001);
Mesh = MeshUpgrade(Mesh,'quadratic');

E = 1; nu = 0; rho = 1; f = {0,0}; gD = {0,0}; gN = {0,0};
function res = u0Func(xy)
  res = 0.1*cos(8*(xy(:,1)+1)).^2.*(8*(xy(:,1)+1)<pi/2);
endfunction
u0 = {'u0Func',0}; v0 = {0,0}; t0 = 0; tend = 4; steps = [4,200];
[u1_all,u2_all,t] = PlaneStressDynamic(Mesh,E,nu,rho,f,gD,gN,u0,v0,t0,tend,steps);

Amp = 0.065; Levels = Amp*[-1:0.1/2:1]; Levels(21)=[]; %% drop Levels = 0
for jj = 2:length(t)
  u1 = u1_all(:,jj); u2 = u2_all(:,jj);
  figure(20+jj); FEMtrimesh(Mesh,u1); zlim(Amp*[-1,1])
                 xlabel('x'); ylabel('y'); zlabel('u_1');
  figure(30+jj); FEMtrimesh(Mesh,u2); zlim(Amp*[-1,1])
                 xlabel('x'); ylabel('y'); zlabel('u_2');
  figure(40+jj); clf; FEMtricontour(Mesh,u1,Levels); axis equal;
                 hold on; plot([x;x(1)],[y;y(1)],'k')
                 xlabel('x'); ylabel('y'); title('u_1');
   figure(50+jj); clf; FEMtricontour(Mesh,u2,Levels); axis equal
                  hold on; plot([x;x(1)],[y;y(1)],'k')
                  xlabel('x'); ylabel('y'); title('u_2')
endfor
