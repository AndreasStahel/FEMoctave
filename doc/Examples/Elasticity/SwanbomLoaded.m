## -*- texinfo -*-
## @deftypefn  {} {} SwanbomLoaded.m
##
## This is a demo file  inside the `doc/Examples/Elasticity/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

%% based on www.youtube.com/watch?v=AXL4ll3aRB8

%% construct the domain
N = 20; R = 0.5e-3; H = 0.5e-3;
al90 = linspace(0,pi/2,N+1)'; xr = 0.5*cos(al90); yr = 0.5*sin(al90);
al360 = linspace(0,2*pi,6*N+1)'; al360 = al360(1:end-1);
Segment.name = 'Segment'; Radius = 0.5;
Segment.border = 1e-3*[2+Radius*cos(al360), Radius*sin(al360), zeros(size(al360))];
x = 1e-3*[0;0;4.5-xr;6;6;4.5-flipud(xr)];
y = 1e-3*[-1.5;1.5;1.5-yr;1;-1;-1.5+flipud(yr)];
BC = -[32;22;22*ones(size(al90));11;22;22*ones(size(al90))];
figure(1); plot(x,y,'-',Segment.border(:,1),Segment.border(:,2))

Mesh = CreateMeshTriangle("Swanbom",[x,y,BC],0.1e-6,Segment);
figure(2); FEMtrimesh(Mesh); axis equal
%%Mesh = MeshUpgrade(Mesh,'quadratic');
Mesh = MeshUpgrade(Mesh,'cubic');

Load = -10;  BoundaryLoad = 0*Load/(3e-3*H);  DiscLoad = 1*Load/(pi*R^2*H);
DeltaT = 0;  %% use 0 fir no thermal prestress or 10 for thermal prestress
E_func = @(xy)E*(1 - 0.25*(((xy(:,1)-2e-3).^2+xy(:,2).^2)<R^2));
fx_func     = @(xy)(DiscLoad*(((xy(:,1)-2e-3).^2+xy(:,2).^2)<R^2));
alphaDeltaT = @(xy)(alpha*DeltaT*(((xy(:,1)-2e-3).^2+xy(:,2).^2)<R^2));

%% solve the plane stress problem
[u1,u2] = PlaneStress(Mesh,E_func,nu,{fx_func,0},{0,0},{BoundaryLoad,0},'thermal',alphaDeltaT);
figure(11); FEMtrimesh(Mesh,u1); xlabel('x'); ylabel('y'); zlabel('u_x'); view([-30,30])
figure(12); FEMtrimesh(Mesh,u2); xlabel('x'); ylabel('y'); zlabel('u_y'); view([-130,45])

[sigma_x,sigma_y,tau_xy] = EvaluateStress(Mesh,u1,u2,E,nu,'thermal',alphaDeltaT);
[eps_xx,eps_yy,eps_xy]   = EvaluateStrain(Mesh,u1,u2);
figure(13); FEMtrimesh(Mesh,eps_xx); xlabel('x'); ylabel('y'); zlabel('\epsilon_{xx}')
figure(14); FEMtrimesh(Mesh,eps_yy); xlabel('x'); ylabel('y'); zlabel('\epsilon_{yy}')

vonMises = EvaluateVonMises(sigma_x,sigma_y,tau_xy)*1e-6;
figure(21); FEMtrimesh(Mesh,vonMises); xlabel('x'); ylabel('y'); zlabel('von Mises [MPa]'); view([35,35])
figure(22); FEMtricontour(Mesh,vonMises,51); xlabel('x'); ylabel('y'); title('von Mises [MPa]')
            colorbar(); axis equal
            hold on; plot([x;x(1)],[y;y(1)],'k'); hold off

figure(23); FEMtrimesh(Mesh,sigma_x*1e-6); xlabel('x'); ylabel('y'); zlabel('\sigma_x [MPa]'); view([35,35])
figure(25); FEMtricontour(Mesh,sigma_x*1e-6,51); xlabel('x'); ylabel('y'); title('\sigma_x [MPa]');
            axis equal; colorbar(); hold on; plot([x;x(1)],[y;y(1)],'k'); hold off
figure(24); FEMtrimesh(Mesh,sigma_y*1e-6); xlabel('x'); ylabel('y'); zlabel('\sigma_y [MPa]')
Max_vonMises = max(vonMises)

phi = linspace(0,2*pi,2001)'; Rc = 0.5001e-3; phiD = phi/pi*180;
x = 0.002+Rc*cos(phi); y = Rc*sin(phi);
[sigma_xc,sigma_yc,tau_xyc] = EvaluateStress(Mesh,u1,u2,E,nu,'thermal',alphaDeltaT,'curve',[x,y]);
sigma_n = sigma_xc.*cos(phi).^2 + 2*tau_xyc.*cos(phi).*sin(phi) + sigma_yc.*sin(phi).^2;
figure(31); plot(phiD,sigma_n*1e-6,phiD,sigma_xc*1e-6); xlabel('angle'); ylabel('\sigma [MPa]')
            legend('\sigma_n','\sigma_x','location','north'); xlim([0,360])
