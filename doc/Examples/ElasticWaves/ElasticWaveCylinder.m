## -*- texinfo -*-
## @deftypefn  {} {} ElasticWaveCylinder.m
##
## This is a demo file  inside the `doc/Examples/ElasticWaves/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

E = 1; nu = 0; rho = 1; f = {0,0}; gD = {0,0}; gN = {0,0}; L = 10; H = L;
function res = u0Func(xy)
  r = sqrt(xy(:,1).^2+xy(:,2).^2);
  res = 0.1*cos(4*r).^2.*(4*r<pi/2);
endfunction
u0 = {'u0Func',0}; v0 = {0,0}; t0 = 0; tend = 4; steps = [4,100];

N = 100;
Mesh = CreateMeshRect(linspace(-L/2,+L/2,N+1),linspace(-H/2,+H/2,N+1),-11,-11,-11,-11);
Mesh = MeshUpgrade(Mesh,'quadratic');
[u1_all,u2_all,t] = PlaneStressDynamic(Mesh,E,nu,rho,f,gD,gN,u0,v0,t0,tend,steps);

Amp = 0.02; Levels = Amp*[-1:0.1:1]; Levels(11)=[];  %% drop Levels = 0
for jj = 2:length(t)
  u1 = u1_all(:,jj); u2 = u2_all(:,jj);
  figure(20+jj); FEMtrimesh(Mesh,u1); zlim(Amp*[-1,1])
                 xlabel('x'); ylabel('y'); zlabel('u_1');
  figure(30+jj); FEMtrimesh(Mesh,u2); zlim(Amp*[-1,1])
                 xlabel('x'); ylabel('y'); zlabel('u_2');
  figure(40+jj); clf; FEMtricontour(Mesh,u1,Levels); axis equal;
                 xlabel('x'); ylabel('y'); title('u_1');
  figure(50+jj); clf; FEMtricontour(Mesh,u2,Levels); axis equal
                 xlabel('x'); ylabel('y'); title('u_2')
endfor


