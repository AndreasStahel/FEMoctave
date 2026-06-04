## -*- texinfo -*-
## @deftypefn  {} {} HeliumGroundState.m
##
## This is a demo file inside the `doc/Examples/Schroedinger/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

pkg load femoctave
R = 100; N = 40;
grid = exp(4*linspace(0,1,N))-1;
grid = R*grid/max(grid);

FEMmesh = CreateMeshRect(grid,grid,-1,-1,-1,-1);
if 0  FEMmesh = MeshUpgrade(FEMmesh,'quadratic');
else  FEMmesh = MeshUpgrade(FEMmesh,'cubic');
endif

a = 1;
b = @(xy) (1*min(xy,[],2)-2*sum(xy,2))./(xy(:,1).*xy(:,2));
%% eigenvalues closest to -2, i.e. most negative
[Eval,Evec] = BVP2Deig(FEMmesh,a,b(FEMmesh.GP),1,0,2,'mode',-2);
Eval

figure(1); FEMtrimesh(FEMmesh,Evec(:,1))
           xlabel('\rho_1'); ylabel('\rho_2')
           xlim([0,8]); ylim([0,8])
figure(2); FEMtricontour(FEMmesh,Evec(:,1),21)
           xlabel('\rho_1'); ylabel('\rho_2')
           xlim([0,8]); ylim([0,8]); axis equal

