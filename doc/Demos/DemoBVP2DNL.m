## -*- texinfo -*-
## @deftypefn  {} {} DemoBVP2DNL.m
##
## This is a demo file inside the `doc/Demos/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

Mesh = CreateMeshTriangle('triangle',[-1,0,-1;1,0,-2;0.5,1,-1],0.01);
Mesh = MeshUpgrade(Mesh,'quadratic');
a = 1; b0 = 0; bx = 0; by = 0; gD = 1; u0 = 1; c = 1;
f = {@(xy,u)c*sin(u);   @(xy,u)c*cos(u)};

u = BVP2DNL(Mesh,a,b0,bx,by,f,gD,1,-1,u0,'tol',1e-10,'Display','iter');
figure(1); FEMtrimesh(Mesh,u);
           xlabel('x'); ylabel('y'); zlabel('u'); view([-80,50])
figure(2); FEMtricontour(Mesh,u); xlabel('x'); ylabel('y'); colorbar
