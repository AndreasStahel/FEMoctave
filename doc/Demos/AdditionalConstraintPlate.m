## -*- texinfo -*-
## @deftypefn  {} {} AdditionalConstraintPlate.m
##
## This is a demo file  inside the `doc/Demos/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

Mesh = CreateMeshRect(linspace(-1,1,9),linspace(-1,1,9),-23, -23, -32, -32);
Mesh = MeshUpgrade(Mesh,'quadratic');
function res = fx(xy)
  res = 1*xy(:,1).*cos(xy(:,2));
endfunction
function res = fy(xy)
  res = 1*xy(:,2).*cos(xy(:,1));
endfunction
[u1,u2] = PlaneStress(Mesh,1,0,{0,0},{0,0},{'fx','fy'});
figure(2); FEMtrimesh(Mesh,u1); xlabel('x'); ylabel('y'); zlabel('u_1')
figure(3); FEMtrimesh(Mesh,u2); xlabel('x'); ylabel('y'); zlabel('u_2')
figure(4); FEMtrimesh(Mesh,sqrt(u1.^2+u2.^2)); xlabel('x'); ylabel('y'); zlabel('|u|')
figure(5); clf; FEMtricontour(Mesh,sqrt(u1.^2+u2.^2)); xlabel('x'); ylabel('y');  axis equal

Pos = [0.0]; Mode = [-1,-1];               %% fix the origin
Mesh = MeshAddConstraint(Mesh,Pos,Mode);   %% remove rotations
Pos = [1,0]; Mode = [-2,-1];
Mesh = MeshAddConstraint(Mesh,Pos,Mode);
[u1m,u2m] = PlaneStress(Mesh,1,0,{0,0},{0,0},{'fx','fy'});

figure(12); FEMtrimesh(Mesh,u1m); xlabel('x'); ylabel('y'); zlabel('u_1')
figure(13); FEMtrimesh(Mesh,u2m); xlabel('x'); ylabel('y'); zlabel('u_2')
figure(14); FEMtrimesh(Mesh,sqrt(u1m.^2+u2m.^2)); xlabel('x'); ylabel('y'); zlabel('|u|')
figure(15); clf; FEMtricontour(Mesh,sqrt(u1m.^2+u2m.^2)); xlabel('x'); ylabel('y'); axis equal

printing = 0;
if printing
  figure(5); print -dpdfcrop AdditionalConstraintPlate1.pdf
  figure(15); print -dpdfcrop AdditionalConstraintPlate2.pdf
endif


