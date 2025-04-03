## -*- texinfo -*-
## @deftypefn  {} {} Liouville.m
##
## This is a demo file inside the `doc/Demos/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

R = 1; N = 36; alpha = linspace(0,pi,N)';
xy = [R*cos(alpha),R*sin(alpha),[-ones(N-1,1);-2]];
Mesh = CreateMeshTriangle('circle',xy,0.01);
%Mesh = MeshUpgrade(Mesh,'quadratic');
Mesh = MeshUpgrade(Mesh,'cubic');

a = 1; b0 = 0; bx = 0; by = 0; u0 = 0; c = 1.0; gD = +0.2;
f = {@(xy,u)exp(c*u);   @(xy,u)c*exp(c*u)};
disp('solving the 2D problem')
u = BVP2DNL(Mesh,a,b0,bx,by,f,gD,0,0,u0,'tol',1e-10,'Display','iter');
figure(1); FEMtrimesh(Mesh,u); xlabel('x'); ylabel('y'); zlabel('u')

Interval = linspace(0,R,201);
f_r = {@(r,u)r.*exp(c*u), @(r,u)c*r.*exp(c*u)};
disp('solving the 1D problem')
[r,u1] = BVP1DNL(Interval,@(r)r,0,0,1,f_r,[0,0],gD,u0,'tol',1e-10,'Display','iter');
u2d = FEMgriddata(Mesh,u,r,zeros(size(r)));
figure(2); plot(r,u1,r,u2d); xlabel('r'); ylabel('u')
           legend('1D radial','2D slice')
figure(3); plot(r,u1-u2d); xlabel('r'); ylabel('difference of u')

