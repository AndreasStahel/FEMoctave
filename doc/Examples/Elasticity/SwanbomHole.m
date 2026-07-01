## -*- texinfo -*-
## @deftypefn  {} {} SwanbomHole.m
##
## This is a demo file  inside the `doc/Examples/Elasticity/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

%% based on www.youtube.com/watch?v=AXL4ll3aRB8

%% construct the domain
N = 10;
al90 = linspace(0,pi/2,N+1)'; xr = 0.5*cos(al90); yr = 0.5*sin(al90);
al360 = linspace(0,2*pi,4*N+1)'; al360 = al360(1:end-1);
Hole.name = 'Hole'; Hole.point = 1e-3*[2,0]; Radius = 0.5e-3;
Hole.border = [2e-3+Radius*cos(al360), Radius*sin(al360), -22*ones(size(al360))];
x = 1e-3*[0;0;4.5-xr;6;6;4.5-flipud(xr)];
y = 1e-3*[-1.5;1.5;1.5-yr;1;-1;-1.5+flipud(yr)];
BC = [32;22;22*ones(size(al90));11;22;22*ones(size(al90))];
figure(1); plot(x,y,'-',Hole.border(:,1),Hole.border(:,2)); axis equal

Mesh = CreateMeshTriangle("Swanbom",[x,y,-BC],0.1e-6,Hole);
figure(2); FEMtrimesh(Mesh); axis equal
%Mesh = MeshUpgrade(Mesh,'quadratic');
Mesh = MeshUpgrade(Mesh,'cubic');

Load = 10;
E = 200e9; nu = 0.25; gN = {-Load/(3*0.5)*1e6,0}; %% data for steel
[u1,u2] = PlaneStress(Mesh,E,nu,{0,0},{0,0},gN);  %% solve the plane stress problem

figure(11); FEMtrimesh(Mesh,u1); xlabel('x'); ylabel('y'); zlabel('u_x'); view([-30,30])
figure(12); FEMtrimesh(Mesh,u2); xlabel('x'); ylabel('y'); zlabel('u_y'); view([-130,45])

[eps_xx,eps_yy,eps_xy]   = EvaluateStrain(Mesh,u1,u2);
figure(13); FEMtrimesh(Mesh,eps_xx); xlabel('x'); ylabel('y'); zlabel('\epsilon_{xx}')
figure(14); FEMtrimesh(Mesh,eps_yy); xlabel('x'); ylabel('y'); zlabel('\epsilon_{yy}')

[sigma_x,sigma_y,tau_xy] = EvaluateStress(Mesh,u1,u2,E,nu);
vonMises = EvaluateVonMises(sigma_x,sigma_y,tau_xy)*1e-6;
figure(21); FEMtrisurf(Mesh,vonMises); xlabel('x'); ylabel('y'); zlabel('von Mises [MPa]')
figure(22); FEMtricontour(Mesh,vonMises,51); xlabel('x'); ylabel('y'); title('von Mises [MPa]')
            colorbar(); axis equal
            hold on; plot([x;x(1)],[y;y(1)],'k'); hold off
figure(23); FEMtrimesh(Mesh,sigma_x*1e-6); xlabel('x'); ylabel('y'); zlabel('\sigma_x [MPa]')
figure(24); FEMtrimesh(Mesh,sigma_y*1e-6); xlabel('x'); ylabel('y'); zlabel('\sigma_y [MPa]')

Max_vonMises = max(vonMises)
