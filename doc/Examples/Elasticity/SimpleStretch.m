## -*- texinfo -*-
## @deftypefn  {} {} SimpleStretch.m
##
## This is a demo file inside the `doc/Examples/Elasticity/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn


if 0  %% only one point constraint, no Saint Venant effect if nu>0
  Mesh = CreateMeshRect(linspace(0,1,21),linspace(0,0.1,17),-22,-22,-12,-22);
  Mesh = MeshAddConstraint(Mesh,[0,0.05],[-1,-1]);
else  %% edge at x=0 constraint, Saint Venant if nu>0
  Mesh = CreateMeshRect(linspace(0,1,21),linspace(0,0.1,11),-22,-22,-11,-22);
endif

%%Mesh = MeshUpgrade(Mesh,'quadratic');
%%Mesh = MeshUpgrade(Mesh,'cubic');
E = 1; nu = 0*0.2;
[u1,u2]= PlaneStress(Mesh,E,nu,{0.01,0.0},{0,0},{0,0});

figure(1); FEMtrimesh(Mesh,u1); xlabel('x'); ylabel('y'); zlabel('u_1')
figure(2); FEMtrimesh(Mesh,u2); xlabel('x'); ylabel('y'); zlabel('u_2')
MaxAbs_u1_u2 = [max(abs(u1)), max(abs(u2))]

[eps_xx,eps_yy,eps_xy] = EvaluateStrain(Mesh,u1,u2);
figure(3); FEMtrimesh(Mesh,eps_xx); xlabel('x'); ylabel('y'); zlabel('\epsilon_{xx}')
figure(4); FEMtrimesh(Mesh,eps_yy); xlabel('x'); ylabel('y'); zlabel('\epsilon_{yy}')
figure(5); FEMtrimesh(Mesh,eps_xy); xlabel('x'); ylabel('y'); zlabel('\epsilon_{xy}')

y = linspace(0,0.1,1001)'; x = 1-2.5e-2*ones(size(y));
[eps_xx_c,eps_yy_c,eps_xy_c] = EvaluateStrain(Mesh,u1,u2,[x,y]);
figure(13); plot(y,eps_xx_c); xlabel('y'); ylabel('\epsilon_{xx}')
figure(14); plot(y,eps_yy_c); xlabel('y'); ylabel('\epsilon_{yy}')
figure(15); plot(y,eps_xy_c); xlabel('y'); ylabel('\epsilon_{xy}')
%% analyse the numerical problem
[~,~,eps_xy1] = FEMgriddata(Mesh,u1,x,y);
[~,eps_xy2,~] = FEMgriddata(Mesh,u2,x,y);
figure(16); plot(y,eps_xy1,y,-eps_xy2);xlabel('y'); %%ylabel('contributions to \epsilon_{xy}')
            legend('+du1/dy','-du2/dx','location','north')

