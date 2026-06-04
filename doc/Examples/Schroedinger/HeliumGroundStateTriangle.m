## -*- texinfo -*-
## @deftypefn  {} {} HeliumGroundStateTriangle.m
##
## This is a demo file inside the `doc/Examples/Schroedinger/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

N = 40;  %% select size of grid

function xy_new = Deform(xy)
  xy_new = exp(4*xy)-1;
  xy_new = xy_new/max(xy_new(:))*100;  %% R = 100
endfunction
FEMmesh = CreateMeshTriangle('Tria',[0,0,-1;1,0,-1;1,1,-2],0.4*(1/N)^2);
FEMmesh = MeshDeform(FEMmesh,'Deform');


if 0  FEMmesh = MeshUpgrade(FEMmesh,'quadratic');
else  FEMmesh = MeshUpgrade(FEMmesh,'cubic');
endif

a = 1; b = @(xy) -(2*xy(:,1)+xy(:,2))./(xy(:,1).*xy(:,2));
[Eval,Evec] = BVP2Deig(FEMmesh,a,b(FEMmesh.GP),1,0,1,'mode',-2);
Eval

[Max,MaxInd] = max(abs(Evec));  %% assure positive maximal value
Evec = sign(Evec(MaxInd))*Evec;

figure(101); FEMtrimesh(FEMmesh, Evec);       xlabel('\rho_1'); ylabel('\rho_2')
             xlim([0,8]); ylim([0,8])
figure(102); FEMtricontour(FEMmesh, Evec,21); xlabel('\rho_1'); ylabel('\rho_2')
             xlim([0,8]); ylim([0,8]); axis equal
