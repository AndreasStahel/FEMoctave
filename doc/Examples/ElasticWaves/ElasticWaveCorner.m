%%% speed sqrt(E/rho) for longitudinal waves
E = 1; nu = 0; rho = 1;
f = {0,0}; gD = {0,0}; gN = {0,0};
function res = u0Func(xy)
  res = 0.1*cos(4*xy(:,1)).^2.*(4*xy(:,1)<pi/2);
endfunction
u0 = {'u0Func',0}; v0 = {0,0};
t0 = 0; tend = 10; steps = [5,200];
L = 12; H = 1;

Mesh = CreateMeshRect(linspace(0,L,121),linspace(-H/2,+H/2,11),-22,-22,-22,-5);
function res = Deform(xy)
  alpha = pi/8;  R = 5/alpha;
  x = xy(:,1); y = xy(:,2) + R; angles = x/R; r = y;
  res = [r.*sin(angles),R-r.*cos(angles)];
endfunction
Mesh = MeshDeform(Mesh,'Deform');
Mesh = MeshUpgrade(Mesh,'quadratic');
solver = 'implicit';
[u1_all,u2_all,t] = PlaneStressDynamic(Mesh,E,nu,rho,f,gD,gN,u0,v0,t0,tend,steps,'solver',solver);

Amp = 0.07;
for jj = 2:length(t)
  u1 = u1_all(:,jj); u2 = u2_all(:,jj);
  disp(sprintf('at time t=%i, max(u1) = %g, max(u2) = %g, max(u) = %g',...
                  t(jj),max(u1),max(u2),max(sqrt(u1.^2+u2.^2))))
  figure(20+jj); FEMtrimesh(Mesh,u1); zlim(Amp*[-0.1,1])
                 xlabel('x'); ylabel('y'); zlabel('u_1');
  figure(30+jj); FEMtrimesh(Mesh,u2); zlim(Amp*[-0.1,1])
                 xlabel('x'); ylabel('y'); zlabel('u_2');
  figure(40+jj); FEMtrimesh(Mesh,sqrt(u1.^2+u2.^2)); zlim(Amp*[-0.1,1])
                 xlabel('x'); ylabel('y'); zlabel('|u|');
endfor


