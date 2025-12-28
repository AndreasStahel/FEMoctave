%% test the IBVP2DNL code
clear *
N = 25;  L = 50;
FEMmesh = CreateMeshRect(linspace(0,L,N),linspace(0,L,N),-2,-1,-2,-1);
FEMmesh = MeshUpgrade(FEMmesh,'quadratic');

m = 1; a = 1; b0 = 0; bx = 0; by = 0;
%% a = [1 1/2 0];  %% for nonistropic diffusion
gD = 0; gN1 = 0; gN2 = 0;
t0 = 0; tend = 20; steps = [25,10];
f = @(xy,t,u)u.*(1-u);
u0 = 0.02*exp(-(FEMmesh.nodes(:,1).^2+FEMmesh.nodes(:,2).^2));

tic()
[u_dyn,t] = IBVP2DNL(FEMmesh,m,a,b0,bx,by,f,gD,gN1,gN2,u0,t0,tend,steps,'solver','CNPC');
toc()

u_end = u_dyn(:,end);
figure(1); FEMtrimesh(FEMmesh,u_end);    xlabel('x'); ylabel('y'); view([70,30])
figure(2); FEMtricontour(FEMmesh,u_end); xlabel('x'); ylabel('y');
           xlim([0,L]); ylim([0,L]); axis equal

if 1 %% animation
  figure(99)
  for jj = 1: size(u_dyn,2)
    FEMtrimesh(FEMmesh,u_dyn(:,jj));  xlabel('x'); ylabel('y'); view([50,30]); zlim([0,1]);
    drawnow(); pause(0.3)
  endfor
endif


printing = 0;
if printing
  figure(1); print -dpdfcrop Fisher2D.pdf
  figure(2); print -dpdfcrop Fisher2Dcontour.pdf
endif

