## -*- texinfo -*-
## @deftypefn  {} {} LaplaceRectangle.m
##
## This is a demo file inside the `doc/Demos/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

FEMmesh = CreateMeshRect([0:0.1:1],[0:0.1:2],-1,-2,-1,-2);
FEMmesh = MeshUpgrade(FEMmesh,'cubic');  %% uncomment to use quadratic elements

u = BVP2Dsym(FEMmesh,1,0,0.25,0,0,0);

figure(1); FEMtrimesh(FEMmesh,u);
           xlabel('x'); ylabel('y');
