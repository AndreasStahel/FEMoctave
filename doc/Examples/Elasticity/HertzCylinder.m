PHASE = 1  %% 1: setup and computation
           %% 2: visualization of the result
           %% 3: parametric study for small penetration depth D
           %% 4: parametric study for large penetration depth D
switch PHASE
case 1  %% PHASE 1
clear *
pkg load femoctave
E = 200e3; nu = 0.24 ;%% N/mm^2  parameter for steel
global R D
R = 1; D = 0.01;          %% radius of cylinder and indentation depth
W = 0.75; H = 1;          %% width and height of the computational domain
a0 = sqrt(D)*sqrt(2*R-D)  %% first estimate of contact point

Area = 0.1^2; MeshType = 'quadratic';  %% definition of the mesh
Seg1.name = 'Segment';
Seg1.border = [0,-2*a0,0;2*a0,-2*a0,0;2*a0,0,0];
Point1.name = 'MeshSize';   Point1.where = [+a0,-a0];  Point1.area = Area/20;   %% finer mesh arounf origin
Point2.name = 'MeshSize'; Point2.where = [+2*a0,-2*a0]; Point2.area = Area;

dd = 0.001;   %% very fine mesh at contact point
Mesh = CreateMeshTriangle('Flat',[0,0,-21;a0-dd,0,-21;a0,0,-22;a0+dd,0,-22;W,0,-12;W,-H,-11;0,-H,-12],...
          Area,Seg1,Point1,Point2);

function res = disp_y(xy)   %% displacement in y-direction as function of x
  global R D
  a0 = sqrt(D)*sqrt(2*R-D);
  res = (R-D-sqrt(R^2 - xy(:,1).^2));
  res = res.*(xy(:,1)<=a0).*(xy(:,2)>-eps);
endfunction

function res = FindFirstPositive(x,sigma_y);  %% find the first zero of the normal stress
  ind = find(sigma_y>0,1);
  x0 = x(ind-1); dx = x(ind)-x0;
  y0 = sigma_y(ind-1); y1 = sigma_y(ind);
  res = x0 -y0*dx/(y1-y0);
endfunction

figure(1); FEMtrimesh(Mesh); axis equal; xlabel('x'); ylabel('y')
           xlim(0.145+[-0.03,0.03]);ylim([-0.05,0.01])
Mesh = MeshUpgrade(Mesh,MeshType);
[u1,u2] = PlaneStrain(Mesh,E,nu,{0,0},{0,'disp_y'},{0,0});

[sigma_x,sigma_y,tau_xy,sigma_z] = EvaluateStress(Mesh,u1,u2,E,nu);
x_edge = linspace(0,2*a0,1000)'; sigma_y_edge = FEMgriddata(Mesh,sigma_y,x_edge,0*x_edge);
figure(30); plot(x_edge,sigma_y_edge,a0,0,'+'); xlabel('x'); ylabel('\sigma_y'); xlim([0,max(x_edge)])
a = FindFirstPositive(x_edge,sigma_y_edge)

for jj = 1:5   %% use 5 iterations and observe the value of a
  Mesh = CreateMeshTriangle('Flat',
         [0,0,-21;a-dd,0,-21;a,0,-22;a+dd,0,-22;W,0,-12;W,-H,-11;0,-H,-12],...
          Area,Seg1,Point1,Point2);
  Mesh = MeshUpgrade(Mesh,MeshType);
  [u1,u2] = PlaneStrain(Mesh,E,nu,{0,0},{0,'disp_y'},{0,0});
  [sigma_x,sigma_y,tau_xy,sigma_z] = EvaluateStress(Mesh,u1,u2,E,nu);
  sigma_y_edge = FEMgriddata(Mesh,sigma_y,x_edge,0*x_edge);
  a = FindFirstPositive(x_edge,sigma_y_edge)
endfor

%% display the normal stress as function of x along the upper edge
ind = find(x_edge<a); sigma_y_fit = mean(sigma_y_edge(ind))*4/pi*sqrt(1-x_edge(ind).^2/a^2);
figure(31); plot(x_edge,sigma_y_edge,x_edge(ind),sigma_y_fit,a0,0,'+',a,0,'+');
            xlabel('x'); ylabel('\sigma_y'); legend('FE data','theory','location','southeast'); xlim([0,1.1*a0])
VonMises = EvaluateVonMises(sigma_x,sigma_y,tau_xy,sigma_z);

%% display the displacement of the upper edge
y_edge = FEMgriddata(Mesh,u2,x_edge,0*x_edge);
y_circle = R-D-sqrt(R^2-x_edge.^2);
figure(20); plot(x_edge,y_edge,x_edge,y_circle,a,0,'*'); xlabel('x'); ylabel('u_2');
            legend('solid','circle', 'location','northwest'); xlim([0,max(x_edge)])

x_edge_low = linspace(0,W,1000)';
sigma_y_edge_low  = FEMgriddata(Mesh,sigma_y,x_edge_low,-H*ones(size(x_edge_low)));
sigma_y_edge_top  = FEMgriddata(Mesh,sigma_y,x_edge_low,zeros(size(x_edge_low)));
sigma_y_edge_half = FEMgriddata(Mesh,sigma_y,x_edge_low,-H/2*ones(size(x_edge_low)));
sigma_y_edge_10 = FEMgriddata(Mesh,sigma_y,x_edge_low,-H/10*ones(size(x_edge_low)));
sigma_y_edge_100 = FEMgriddata(Mesh,sigma_y,x_edge_low,-H/100*ones(size(x_edge_low)));
figure(90); plot(x_edge_low,sigma_y_edge_top,x_edge_low,sigma_y_edge_100,x_edge_low,sigma_y_edge_10,...
            x_edge_low,sigma_y_edge_half,x_edge_low,sigma_y_edge_low);
            xlabel('x'); ylabel('\sigma_y'); xlim([0,max(x_edge)])
            legend('at y = 0','at y = -H/100','at y = -H/10','at y = -H/2','at y = -H','location','southeast')

PressureUpperEdgeLocal = trapz(x_edge,sigma_y_edge)
PressureUpperEdge = trapz(x_edge_low,sigma_y_edge_top)
Pressure100Edge = trapz(x_edge_low,sigma_y_edge_100)
Pressure10Edge = trapz(x_edge_low,sigma_y_edge_10)
PressureHalfEdge  = trapz(x_edge_low,sigma_y_edge_half)
PressureLowerEdge = trapz(x_edge_low,sigma_y_edge_low)

Estar = E/(1-nu^2);
P = -2*PressureLowerEdge;
a = sqrt(4*P*R/(pi*Estar))

case 2  %% PHASE 2

[x,y] = meshgrid(linspace(0,3*a,51),linspace(-3*a,0,51));
u1g = FEMgriddata(Mesh,u1,x,y); u2g = FEMgriddata(Mesh,u2,x,y);
x_g = x + u1g;  y_g = y + u2g;              %% construct the deformed grid
sigma_yg = FEMgriddata(Mesh,sigma_y,x,y);   %% evaluate sigma_y on the  grid
figure(101); mesh(x_g,y_g,sigma_yg); xlabel('x'); ylabel('y');zlabel('\sigma_y')
figure(111); contourf(x_g,y_g,sigma_yg/1e3,linspace(min(sigma_yg(:)),0,51)/1e3); xlabel('x'); ylabel('y');
             title('\sigma_y [kPa]'); colorbar
figure(121); contourf(x_g,y_g,-u2g*1e3,51); xlabel('x'); ylabel('y');
             title('-u_y [\mum]'); colorbar
sigma_xg = FEMgriddata(Mesh,sigma_x,x,y);
figure(102); mesh(x_g,y_g,sigma_xg); xlabel('x'); ylabel('y');zlabel('\sigma_x')
figure(112); contourf(x_g,y_g,sigma_xg/1e3,linspace(min(sigma_yg(:)),0,51)/1e3); xlabel('x'); ylabel('y');
             title('\sigma_x [kPa]'); colorbar
figure(122); contourf(x_g,y_g,u1g*1e3,51); xlabel('x'); ylabel('y');
             title('u_x [\mum]'); colorbar

VonMises_g = FEMgriddata(Mesh,VonMises,x,y);
figure(103); mesh(x,y,VonMises_g); xlabel('x'); ylabel('y');zlabel('von Mises')
figure(113); contourf(x_g,y_g,VonMises_g/1e3,linspace(0,max(VonMises_g(:)),51)/1e3); xlabel('x'); ylabel('y');
             title('von Mises [kPa]'); colorbar
figure(123); contourf(x_g,y_g,sqrt(u1g.^2+u2g.^2)*1e3,51); xlabel('x'); ylabel('y');
             title('displacement [\mum]'); colorbar

DOriginal = D; aOriginal = a;

case 3  %% Phase 3: parametric study for small D
a = aOriginal;
D_List = DOriginal*[1:-0.1:0.1]';
a_List = zeros(size(D_List)); Pressure_List = a_List;
a_List(1) = aOriginal; Pressure_List(1) = -2*PressureLowerEdge;
for ii = 1:length(D_List)
  D = D_List(ii);
  disp(sprintf('working with penetration depth D = %g',D))
  a_old = 0; jj = 0;
  do
    jj++;
    Mesh = CreateMeshTriangle('Flat',
           [0,0,-21;a-dd,0,-21;a,0,-22;a+dd,0,-22;W,0,-12;W,-H,-11;0,-H,-12],...
            Area,Seg1,Point1,Point2);
    Mesh = MeshUpgrade(Mesh,MeshType);
    [u1,u2] = PlaneStrain(Mesh,E,nu,{0,0},{0,'disp_y'},{0,0});
    [sigma_x,sigma_y,tau_xy,sigma_z] = EvaluateStress(Mesh,u1,u2,E,nu);
    sigma_y_edge = FEMgriddata(Mesh,sigma_y,x_edge,0*x_edge);
    a_old = a;
    a = FindFirstPositive(x_edge,sigma_y_edge);
    disp(sprintf('value a = %g, last change = %g', a , a-a_old))
  until or(abs(a-a_old)< 5e-5,jj>20)
  a_List(ii) = a;
  sigma_y_edge_low  = FEMgriddata(Mesh,sigma_y,x_edge_low,-H*ones(size(x_edge_low)));
  PressureLowerEdge = trapz(x_edge_low,sigma_y_edge_low);
  Pressure_List(ii) = -2 * PressureLowerEdge;
endfor
D_List = [D_List;0]; Pressure_List = [Pressure_List;0];  a_List = [a_List;0];
figure(201); plot(Pressure_List,D_List); xlabel('pressure'); ylabel('penetration depth')
figure(202); plot(Pressure_List,a_List); xlabel('pressure'); ylabel('contact width')

p = LinearRegression(a_List.^2,Pressure_List);
figure(203); plot(a_List,Pressure_List,'*',a_List,p*a_List.^2)
             ylabel('pressure'); xlabel('contact width')
             legend('FEM','c*a^2','location','northwest')

case 4  %% Phase 4: parametric study for large D
a = aOriginal;
D_List = DOriginal*[1:0.4:10]';
a_List = zeros(size(D_List)); Pressure_List = a_List;

for ii = 1:length(D_List)
  D = D_List(ii);
  disp(sprintf('working with penetration depth D = %g',D))
  x_edge = linspace(0,1.3*a,1000');
  a_old = 0;jj = 0;
  do
  jj++;
  Mesh = CreateMeshTriangle('Flat',
           [0,0,-21;a-dd,0,-21;a,0,-22;a+dd,0,-22;W,0,-12;W,-H,-11;0,-H,-12],...
            Area,Seg1,Point1,Point2);
    Mesh = MeshUpgrade(Mesh,MeshType);
    [u1,u2] = PlaneStrain(Mesh,E,nu,{0,0},{0,'disp_y'},{0,0});
    [sigma_x,sigma_y,tau_xy,sigma_z] = EvaluateStress(Mesh,u1,u2,E,nu);
    sigma_y_edge = FEMgriddata(Mesh,sigma_y,x_edge,0*x_edge);
    a_old = a;
    a = FindFirstPositive(x_edge,sigma_y_edge);
    disp(sprintf('iteration %i, value a = %g, last change = %g', jj, a, a-a_old))
  until or(abs(a-a_old)< 5e-5,jj>20)
  a_List(ii) = a;
  sigma_y_edge_low  = FEMgriddata(Mesh,sigma_y,x_edge_low,-H*ones(size(x_edge_low)));
  PressureLowerEdge = trapz(x_edge_low,sigma_y_edge_low);
  Pressure_List(ii) = -2 * PressureLowerEdge;
endfor

D_List = [0;D_List]; Pressure_List = [0;Pressure_List];  a_List = [0;a_List];
figure(301); plot(Pressure_List,D_List); xlabel('pressure'); ylabel('penetration depth')
figure(302); plot(Pressure_List,a_List); xlabel('pressure'); ylabel('contact width')

p = LinearRegression(a_List.^2,Pressure_List);
figure(303); plot(a_List,Pressure_List,'*',a_List,p*a_List.^2)
             ylabel('pressure'); xlabel('contact width')
             legend('FEM','c*a^2','location','northwest')
endswitch
