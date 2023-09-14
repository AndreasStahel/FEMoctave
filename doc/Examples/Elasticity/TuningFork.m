angle = linspace(0,pi,11); halfangle = linspace(0,pi/2,6);
R1 = 5.5; R2 = 6.25; H1 = 108.5; H2 = H1+15;
dH = 0; %% for different length of the two branches

x = [R1*cos(fliplr(angle)),R1,R1+7,R1+7-R2*sin(halfangle),R1+7-R2,-R1-7+R2,...
     -R1-7+R2*sin(fliplr(halfangle)),-R1-7,-R1]/1000;
y = [-R1*sin(angle),H1+dH,H1+dH,H1-H2-R2+R2*cos(halfangle),-50,-50,...
     H1-H2-R2+R2*cos(fliplr(halfangle)),H1,H1]/1000;
%%figure(1); plot(x,y); axis equal

CASE = 1
IND = find(y==-50/1000,1);
switch CASE
case 1  %% only lower edge fixed
  xy = [x',y',-22*ones(length(x),1)]; xy(IND,3)=-11;
case 2  %% lower edge and sides fixed
  xy = [x',y',-22*ones(length(x),1)]; xy([IND-1:IND+1],3)=-11;
endswitch

FEMmesh = CreateMeshTriangle('Fork',xy,3e-6);
%figure(2); FEMtrimesh(FEMmesh); axis equal
FEMmesh = MeshUpgrade(FEMmesh,'quadratic');

E = 200e9; nu = 0.21; rho = 7.9e3;  %% use SI units
[la,u1all,u2all] = PlaneStressEig(FEMmesh,E,nu,rho,6);
Frequencies = sqrt(la')/(2*pi)

Mode = 2
u1 = u1all(:,Mode); u2 = u2all(:,Mode);
MaxDisp = max(max(abs(u1)),max(abs(u2)));  %% scale the values
u1 = 0.005*u1/MaxDisp; u2 = 0.005*u2/MaxDisp;

figure(11); FEMtrimesh(FEMmesh,u1); xlabel('x'); ylabel('y'); zlabel('u_1')
figure(12); FEMtrimesh(FEMmesh,u2); xlabel('x'); ylabel('y'); zlabel('u_2')

[sigma_x,sigma_y,tau_xy] = EvaluateStress(FEMmesh,u1,u2,E,nu);
figure(21); FEMtrimesh(FEMmesh,sigma_x/1e6); xlabel('x'); ylabel('y'); zlabel('\sigma_x')
figure(22); FEMtrimesh(FEMmesh,sigma_y/1e6); xlabel('x'); ylabel('y'); zlabel('\sigma_y')
figure(23); FEMtrimesh(FEMmesh,tau_xy/1e6); xlabel('x'); ylabel('y'); zlabel('\tau_{xy}')
figure(30); ShowDeformation(FEMmesh,u1,u2,1); axis equal

