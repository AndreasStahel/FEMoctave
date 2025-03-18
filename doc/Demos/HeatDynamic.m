%% a dynamic heat problem
clear *
N = 25; %%grid size 2N+1 by 2N+1
%% generate the mesh
if 0  %% first order elements
  FEMmesh = CreateMeshRect(linspace(0,2,2*N+1),linspace(0,2,2*N+1),-1,-1,-1,-1);
else  %% second order elements
  FEMmesh = CreateMeshRect(linspace(0,2,N+1),linspace(0,2,N+1),-1,-1,-1,-1);
  FEMmesh = MeshUpgrade(FEMmesh);
endif
x = FEMmesh.nodes(:,1);y = FEMmesh.nodes(:,2);
%% setup and solve the initial boundary value problem
m = 1; a = 1; b0 = 0; bx = 10; by = 5; f = 0.1; gD = 0; gN1 = 0; gN2 = 0;
t0 = 0; tend = 0.1 ; steps = [6,10];
u0 = zeros(length(FEMmesh.nodes),1);
u0 = 0.005*(2-x).^2.*x.*y.*(2-y);
tic()
[u_dyn,t] = IBVP2D(FEMmesh,m,a,b0,bx,by,f,gD,gN1,gN2,u0,t0,tend,steps,'solver','RK');
toc()

figure(1); FEMtrimesh(FEMmesh,u_dyn(:,end))
           xlabel('x'); ylabel('y')

%% show the animation on screen
u_max = max(u_dyn(:));
for t_ii = 1:length(t)
  figure(2); FEMtrimesh(FEMmesh,u_dyn(:,t_ii)); xlabel('x'); ylabel('y')
             caxis([0,u_max]); axis([0 2 0 2 0 u_max])
  drawnow();
  figure(3); FEMtricontour(FEMmesh,u_dyn(:,t_ii),linspace(0,0.99*u_max,11))
             xlabel('x'); ylabel('y'); caxis([0,u_max]);
  drawnow();
  pause(0.4)
endfor
