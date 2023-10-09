PHASE = 1  %% 1: setup and computation
           %% 2: visualization of the result
           %% 3: parametric study for small and large penetration depth D
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
Seg1.border = [0,-4*a0,0;2*a0,-2*a0,0;2*a0,0,0];
Point1.name = 'MeshSize';
Point1.where = [+a0,-a0];
Point1.area = Area/20;   %% finer mesh around origin
Point2.name = 'MeshSize';
Point2.where = [+2*a0,-2*a0];
Point2.area = Area;

dd = 0.001;   %% very fine mesh at contact point
Mesh = CreateMeshTriangle('Flat',[0,0,-21;a0-dd,0,-21;a0,0,-22;a0+dd,0,-22;W,0,-12;W,-H,-11;0,-H,-12],...
          Area,Seg1,Point1,Point2);

function res = disp_z(xy)   %% displacement in y-direction as function of x
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
Mesh = MeshUpgrade(Mesh,MeshType);
[ur,uz] = AxiStress(Mesh,E,nu,{0,0},{0,'disp_z'},{0,0});

[sigma_r,sigma_y,sigma_z,tau_xz] = EvaluateStressAxi(Mesh,ur,uz,E,nu);
r_edge = linspace(0,2*a0,1000)'; sigma_z_edge = FEMgriddata(Mesh,sigma_z,r_edge,0*r_edge);
figure(30); plot(r_edge,sigma_z_edge,a0,0,'+'); xlabel('r'); ylabel('\sigma_z'); xlim([0,max(r_edge)])
a = FindFirstPositive(r_edge,sigma_z_edge)
jj = 0;
do
  jj++;
  Mesh = CreateMeshTriangle('Flat',
         [0,0,-21;a-dd,0,-21;a,0,-22;a+dd,0,-22;W,0,-12;W,-H,-11;0,-H,-12],...
          Area,Seg1,Point1,Point2);
  Mesh = MeshUpgrade(Mesh,MeshType);
  [ur,uz] = AxiStress(Mesh,E,nu,{0,0},{0,'disp_z'},{0,0});
  [sigma_r,sigma_y,sigma_z,tau_xz] = EvaluateStressAxi(Mesh,ur,uz,E,nu);
  sigma_z_edge = FEMgriddata(Mesh,sigma_z,r_edge,0*r_edge);
  a_old = a;
  a = FindFirstPositive(r_edge,sigma_z_edge);
  disp(sprintf('iteration %i, value a = %g, last change = %g', jj, a , a-a_old))
%%  figure(30+jj); plot(x_edge,sigma_y_edge,x_c,0,'+',x_c1,0,'+'); xlabel('x'); ylabel('\sigma_y')
until  or(abs(a-a_old)< 5e-5,jj>15)

%% display the normal stress as function of x along the upper edge
ind = find(r_edge<a); sigma_z_fit = mean(sigma_z_edge(ind))*4/pi*sqrt(1-r_edge(ind).^2/a^2);
Estar = E/(1-nu^2);
P_Hertz  = 4*E/((1-nu^2)*3)*sqrt(R)*D^(3/2)
a_Hertz = (3*P_Hertz*R/(4*Estar))^(1/3)
sigma_z_Hertz = -3*P_Hertz/(2*pi*a_Hertz^2)*sqrt(1-r_edge(ind).^2/(a_Hertz^2));
figure(31); plot(r_edge,sigma_z_edge,r_edge(ind),sigma_z_fit,r_edge(ind),sigma_z_Hertz,a0,0,'+',a,0,'+');
            xlabel('r'); ylabel('\sigma_z'); legend('FE data','fit','Hertz','location','southeast'); xlim([0,1.1*a0])
VonMises = EvaluateVonMisesAxi(sigma_r,sigma_y,sigma_z,tau_xz);

%% display the displacement of the upper edge
z_edge = FEMgriddata(Mesh,uz,r_edge,0*r_edge);
z_circle = R-D-sqrt(R^2-r_edge.^2);
figure(20); plot(r_edge,z_edge,r_edge,z_circle,a,0,'*'); xlabel('r'); ylabel('u_z');
            legend('solid','circle', 'location','northwest'); xlim([0,max(r_edge)])

r_edge_low = linspace(0,W,1000)';
sigma_z_edge_low  = FEMgriddata(Mesh,sigma_z,r_edge_low,-H*ones(size(r_edge_low)));
sigma_z_edge_top  = FEMgriddata(Mesh,sigma_z,r_edge_low,zeros(size(r_edge_low)));
sigma_z_edge_half = FEMgriddata(Mesh,sigma_z,r_edge_low,-H/2*ones(size(r_edge_low)));
sigma_z_edge_10 = FEMgriddata(Mesh,sigma_z,r_edge_low,-H/10*ones(size(r_edge_low)));
sigma_z_edge_100 = FEMgriddata(Mesh,sigma_z,r_edge_low,-H/100*ones(size(r_edge_low)));
figure(90); plot(r_edge_low,sigma_z_edge_top,r_edge_low,sigma_z_edge_100,r_edge_low,sigma_z_edge_10,...
            r_edge_low,sigma_z_edge_half,r_edge_low,sigma_z_edge_low);
            xlabel('r'); ylabel('\sigma_z'); xlim([0,max(r_edge)])
            legend('at z = 0','at z = -H/100','at z = -H/10','at z = -H/2','at z = -H','location','southeast')

PressureUpperEdgeLocal = trapz(r_edge,2*pi*r_edge.*sigma_z_edge)
PressureUpperEdge = trapz(r_edge_low,2*pi*r_edge_low.*sigma_z_edge_top)
Pressure100Edge = trapz(r_edge_low,2*pi*r_edge_low.*sigma_z_edge_100)
Pressure10Edge = trapz(r_edge_low,2*pi*r_edge_low.*sigma_z_edge_10)
PressureHalfEdge  = trapz(r_edge_low,2*pi*r_edge_low.*sigma_z_edge_half)
PressureLowerEdge = trapz(r_edge_low,2*pi*r_edge_low.*sigma_z_edge_low)

Estar = E/(1-nu^2);            %% Hertz theory
F = -PressureLowerEdge;
a_Hertz = (3*F*R/(4*Estar))^(1/3)
D_Hertz = (3*F/(4*Estar))^(2/3)/R^(1/3)
P_Hertz = 4*Estar/3*sqrt(R)*D^(3/2)

case 2  %% PHASE 2
[r,z] = meshgrid(linspace(0,3*a,51),linspace(-3*a,0,51));
urg = FEMgriddata(Mesh,ur,r,z); uzg = FEMgriddata(Mesh,uz,r,z);
r_g = r + urg;  z_g = z + uzg;              %% construct the deformed grid
sigma_zg = FEMgriddata(Mesh,sigma_z,r,z);   %% evaluate sigma_z on the  grid
figure(101); mesh(r_g,z_g,sigma_zg); xlabel('r'); ylabel('z');zlabel('\sigma_z')
figure(111); contourf(r_g,z_g,sigma_zg/1000,linspace(min(sigma_zg(:))/1000,0,51)); xlabel('r'); ylabel('z');
             title('\sigma_z [kPa]'); colorbar
figure(121); contourf(r_g,z_g,-uzg*1e3,51); xlabel('r'); ylabel('z');
             title('-u_z [\mum]'); colorbar
sigma_rg = FEMgriddata(Mesh,sigma_r,r,z);
figure(102); mesh(r_g,z_g,sigma_rg); xlabel('r'); ylabel('z');zlabel('\sigma_r')
figure(112); contourf(r_g,z_g,sigma_rg/1000,linspace(min(sigma_rg(:))/1000,0,51)); xlabel('r'); ylabel('z');
             title('\sigma_r [kPa]'); colorbar
figure(122); contourf(r_g,z_g,urg*1e3,51); xlabel('r'); ylabel('z');
             title('u_r [\mum]'); colorbar

VonMises_g = FEMgriddata(Mesh,VonMises,r,z);
figure(103); mesh(r,z,VonMises_g); xlabel('r'); ylabel('z');zlabel('von Mises')
figure(113); contourf(r_g,z_g,VonMises_g/1000,linspace(0,max(VonMises_g(:))/1000,51)); xlabel('r'); ylabel('z');
             title('von Mises [kPa]'); colorbar
figure(123); contourf(r_g,z_g,sqrt(urg.^2+uzg.^2)*1e3,51); xlabel('r'); ylabel('z');
             title('displacement [\mum]'); colorbar

DOriginal = D; aOriginal = a;

case 3  %% Phase 3: parametric study for small and large D
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
    [ur,uz] = AxiStress(Mesh,E,nu,{0,0},{0,'disp_z'},{0,0});
    [sigma_r,sigma_y,sigma_z,tau_xy] = EvaluateStressAxi(Mesh,ur,uz,E,nu);
    sigma_z_edge = FEMgriddata(Mesh,sigma_z,r_edge,0*r_edge);
    a_old = a;
    a = FindFirstPositive(r_edge,sigma_z_edge);
    disp(sprintf('iteration %i, value a = %g, last change = %g', jj, a, a-a_old))
  until or(abs(a-a_old)< 5e-5,jj>20)
  a_List(ii) = a;
  sigma_z_edge_low  = FEMgriddata(Mesh,sigma_z,r_edge_low,-H*ones(size(r_edge_low)));
  PressureLowerEdge = trapz(r_edge_low,2*pi*r_edge_low.*sigma_z_edge_low);
  Pressure_List(ii) = -PressureLowerEdge;
endfor
D_List = [D_List;0]; Pressure_List = [Pressure_List;0];  a_List = [a_List;0];

%% save the data
a_List_small = flipud(a_List); Pressure_List_small = flipud(Pressure_List);
D_List_small = flipud(D_List);

a = aOriginal;
D_List = DOriginal*[1.4:0.4:10]';
a_List = zeros(size(D_List)); Pressure_List = a_List;

for ii = 1:length(D_List)
  D = D_List(ii);
  disp(sprintf('working with penetration depth D = %g',D))
  r_edge = linspace(0,1.3*a,1000)';
  jj = 0;
  do
  jj++;
  Mesh = CreateMeshTriangle('Flat',
           [0,0,-21;a-dd,0,-21;a,0,-22;a+dd,0,-22;W,0,-12;W,-H,-11;0,-H,-12],...
            Area,Seg1,Point1,Point2);
    Mesh = MeshUpgrade(Mesh,MeshType);
    [ur,uz] = AxiStress(Mesh,E,nu,{0,0},{0,'disp_z'},{0,0});
    [sigma_r,sigma_y,sigma_z,tau_xy] = EvaluateStressAxi(Mesh,ur,uz,E,nu);
    sigma_z_edge = FEMgriddata(Mesh,sigma_z,r_edge,0*r_edge);
    a_old = a;
    a = FindFirstPositive(r_edge,sigma_z_edge);
    disp(sprintf('iteration %i, value a = %g, last change = %g', jj, a, a-a_old))
  until or(abs(a-a_old)< 5e-5,jj>24)
  a_List(ii) = a;
  a = 1.05*a;
  sigma_z_edge_low  = FEMgriddata(Mesh,sigma_z,r_edge_low,-H*ones(size(r_edge_low)));
  PressureLowerEdge = trapz(r_edge_low,2*pi*r_edge_low.*sigma_z_edge_low);
  Pressure_List(ii) = -PressureLowerEdge;
endfor

D_List = [D_List_small;D_List];
Pressure_List = [Pressure_List_small;Pressure_List];
a_List = [a_List_small;a_List];
D_Hertz = (3/(4*Estar))^(2/3)/R^(1/3)*Pressure_List.^(2/3);

p1 = sum(D_Hertz.*D_List)/sum(D_Hertz.^2)
figure(201); plot(Pressure_List,D_List,'*',Pressure_List,p1*D_Hertz,Pressure_List,D_Hertz);
             xlabel('pressure'); ylabel('penetration depth')
             legend('FE data','fit','Hertz','location','southeast')
figure(202); plot(Pressure_List,D_List,'*',Pressure_List,p1*D_Hertz,Pressure_List,D_Hertz);
             xlabel('pressure'); ylabel('penetration depth')
             legend('FE data','fit','Hertz','location','southeast')
             xlim([0,800]); ylim([0,0.02])

p2 = LinearRegression(a_List.^3,Pressure_List);
figure(203); plot(a_List,Pressure_List,'*',a_List,p2*a_List.^3)
             ylabel('pressure'); xlabel('contact width')
             legend('FEM','c*a^3','location','northwest')
endswitch
