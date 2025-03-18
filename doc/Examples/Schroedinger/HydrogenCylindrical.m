## -*- texinfo -*-
## @deftypefn  {} {} HydrogenCylindrical.m
##
## This is a demo file inside the `doc/Examples/Schroedinger/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

%% evaluate energies for the hydrogen atom in cylindrical coordinates
global R
R = 35; N = 50; alpha = linspace(-pi/2,pi/2,N)';
Rshow = 15;
MeshContour1 =  [cos(alpha),sin(alpha),-ones(size(alpha))];
MeshContour1(N,3) = -2;
FEMmesh = CreateMeshTriangle('circle',MeshContour1,0.25e-3);  %% good value with R = 35
function xy_new = Deform(xy)
  global R;
  CC = 5;
  r = sqrt(xy(:,1).^2 + xy(:,2).^2);
  xy_new = xy.*exp(CC*r)/exp(CC)*R;
endfunction
FEMmesh = MeshDeform(FEMmesh,'Deform');
FEMmesh = MeshUpgrade(FEMmesh,'cubic');

function res = rho(rho_z)
  res = rho_z(:,1);
endfunction

global Shift
Shift = 20;

function res = b0(rho_z)
  global Shift
  res = -2*rho_z(:,1)./sqrt(rho_z(:,1).^2+rho_z(:,2).^2) + Shift*rho_z(:,1);
endfunction

Nval = 5;
[Eval,Evec] = BVP2Deig(FEMmesh,'rho','b0','rho',0,Nval);
Eigenvalues = Eval-Shift
Ratios = Eigenvalues(1)./Eigenvalues

%% force the maximal absolute value positive
[Max,Ind] = max(abs(Evec)); Sign = zeros(Nval,1);
for ii=1:Nval
  Sign(ii) = sign(Evec(Ind(ii),ii));
endfor
Evec = Evec*diag(Sign);


n = 3;
u = Evec(:,n);
figure(1); FEMtrimesh(FEMmesh,u)
           xlabel('x  [a_0]') ; ylabel('z  [a_0]'); zlabel('u')
           xlim([0,Rshow]); ylim([-Rshow,+Rshow])

figure(2); clf; FEMtricontour(FEMmesh,u)
           xlabel('x  [a_0]') ; ylabel('z  [a_0]'); axis equal
           hold on; plot([MeshContour(:,1);0],[MeshContour(:,2);-R],'k')
           xlim(1.5*[0,Rshow]); ylim(0.7*[-Rshow,+Rshow])

r = sqrt(FEMmesh.nodes(:,1).^2+FEMmesh.nodes(:,2).^2);
figure(3); FEMtrimesh(FEMmesh,r.*u.^2)
           xlabel('x  [a_0]') ; ylabel('z  [a_0]'); zlabel('PDF')
           xlim([0,Rshow]); ylim([-Rshow,+Rshow])

x = linspace(0,R,200);
u_x = FEMgriddata(FEMmesh,u,zeros(size(x)),x);

figure(4); plot(x,u_x)
           xlabel('\rho  [a_0]'); ylabel('u')
