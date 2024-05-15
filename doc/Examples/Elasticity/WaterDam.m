## -*- texinfo -*-
## @deftypefn  {} {} WaterDam.m
##
## This is a demo file  inside the `doc/Examples/Elasticity/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

clear *
global H Hwater BaseLeft
H = 30; Base = 0.7*H; Crest = 0.2*H; Hwater = 0.9*H; BaseLeft = H*0.2;
E = 20e9; nu = 0.2;
CASE = 3
switch CASE
case 1 %% no crack
  xy = [0,0,-11;Base,0,-22;Crest,H,-22;0,H,-33]; x = [xy(:,1);xy(1,1)]; y = [xy(:,2);xy(1,2)];
case 2 %% with crack
  h = 0.1; depth = 1;
  xy = [0,0,-11;Base,0,-22;Crest,H,-22;0,H,-33;0,H/2+h,-22;depth,H/2,-22;0,H/2-h,-33];
  x = [xy(:,1);xy(1,1)]; y = [xy(:,2);xy(1,2)];
case 3  %% with foot
  BaseLeft = H*0.2;
  xy = [-BaseLeft,0,-11;Base,0,-22;Crest,H,-22;0,H,-33;0,BaseLeft,-33];
  x = [xy(:,1);xy(1,1)]; y = [xy(:,2);xy(1,2)];
case 4  %% with slope on both sides
  BaseLeft = H*0.2;
  xy = [-BaseLeft,0,-11;Base,0,-22;Crest,H,-22;0,H,-33];
  x = [xy(:,1);xy(1,1)]; y = [xy(:,2);xy(1,2)];
endswitch

FEMmesh = CreateMeshTriangle('Dam',xy,1);
FEMmesh = MeshUpgrade(FEMmesh,'cubic');
figure(1); FEMtrimesh(FEMmesh); xlabel('x'); ylabel('h')

function res = f_dam(xy,dummy)
   global H BaseLeft
   rho = 2.4e3;
   res = -9.81*rho*(H-xy(:,2));
endfunction

switch CASE  %% different surface pressures
case {1,2}
 function res = px(xy,dummy)
   global Hwater
   res = +9.81e3*(Hwater-xy(:,2)).*(xy(:,1)<eps).*(xy(:,2)<Hwater);
 endfunction
 function res = ph(xy,dummy)
   global Hwater
   res = +0*9.81e3*(Hwater-xy(:,2)).*(xy(:,1)<eps).*(xy(:,2)<Hwater);
 endfunction

case 3
 function res = px(xy,dummy)
   global Hwater BaseLeft
   res = +9.81e3*(Hwater-xy(:,2)).*(xy(:,1)<eps).*...
              ((xy(:,2)<Hwater).*(xy(:,2)>BaseLeft)+1/sqrt(2)*(xy(:,2)<=BaseLeft));
 endfunction
 function res = ph(xy,dummy)
   global Hwater BaseLeft
   res = +9.81e3*(Hwater-xy(:,2)).*(xy(:,1)<eps).*(xy(:,2)<=BaseLeft)/sqrt(2);
 endfunction

case 4
 function res = px(xy,dummy)
   global Hwater BaseLeft H
   alpha = atan(BaseLeft/H);
   res = +9.81e3*(Hwater-xy(:,2)).*(xy(:,1)<eps).*...
              ((xy(:,2)<Hwater).*(xy(:,2)>BaseLeft)+cos(alpha)*(xy(:,2)<=BaseLeft));
 endfunction
 function res = ph(xy,dummy)
   global Hwater BaseLeft H
   alpha = atan(BaseLeft/H);
   res = +9.81e3*(Hwater-xy(:,2)).*(xy(:,1)<eps).*(xy(:,2)<=BaseLeft)*sin(alpha);
 endfunction
endswitch

[u1,u2] = PlaneStrain(FEMmesh,E,nu,{0,'f_dam'},{0,0},{'px','ph'});
[sigma_x,sigma_y,tau_xy,sigma_z] = EvaluateStress(FEMmesh,u1,u2,E,nu);
vonMises = EvaluateVonMises(sigma_x,sigma_y,tau_xy,sigma_z);


figure(2); FEMtrimesh(FEMmesh,u1);  xlabel('x'); ylabel('h'); zlabel('u_x'); view([-120,20])
figure(3); FEMtrimesh(FEMmesh,u2);  xlabel('x'); ylabel('h'); zlabel('u_h'); view([-120,20])
figure(11); FEMtrimesh(FEMmesh,sigma_y*1e-6); xlabel('x'); ylabel('h'); zlabel('\sigma_h [MPa]'); view([60,20])
figure(21); clf; FEMtricontour(FEMmesh,sigma_y*1e-6); xlabel('x'); ylabel('h'); colorbar();
hold on; plot(x,y,'k'); title('\sigma_h [MPa]')
figure(12); FEMtrimesh(FEMmesh,sigma_x*1e-6); xlabel('x'); ylabel('h'); zlabel('\sigma_x [MPa]'); view([60,20])
figure(22); clf; FEMtricontour(FEMmesh,sigma_x*1e-6); xlabel('x'); ylabel('h'); colorbar()
hold on; plot(x,y,'k'); title('\sigma_x [MPa]')
figure(13); FEMtrimesh(FEMmesh,vonMises*1e-6); xlabel('x'); ylabel('h'); zlabel('von Mises [MPa]'); view([60,20])
figure(23); clf; FEMtricontour(FEMmesh,vonMises*1e-6); xlabel('x'); ylabel('h'); colorbar()
hold on; plot(x,y,'k'); title('von Mises [MPa]')

Max_vonMises = max(vonMises)/1e6
