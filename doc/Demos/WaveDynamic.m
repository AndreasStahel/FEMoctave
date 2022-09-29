%% a dynamic wave problem
clear *
%% generate a circle
alpha = linspace(0,2*pi,101)'; alpha = alpha(1:end-1); R = 6;
xy = [R*cos(alpha),R*sin(alpha),-ones(size(alpha))];
if 1  %% linear elements
  FEMmesh = CreateMeshTriangle('Circle',xy,0.03);
else  %% quadratic elements
  FEMmesh = CreateMeshTriangle('Circle',xy,4*0.03);
  FEMmesh = MeshUpgrade(FEMmesh,'quadratic');
endif

x = FEMmesh.nodes(:,1); y = FEMmesh.nodes(:,2);
v0 = zeros(size(x)); 
u0 = 0.1*exp(-1*((x-1).^2+y.^2)); u0 = u0.*(R^2-x.^2-y.^2)/R^2;
%% setup and solve the initial boundary value problem
m=1; d=0; a=1; b0=0; bx=0; by=0; f=0; gD=0; gN1=0; gN2=0;
t0=0; tend=7 ; steps=[14,30];
tic();
[u_dyn,t] = I2BVP2D(FEMmesh,m,d,a,b0,bx,by,f,gD,gN1,gN2,u0,v0,t0,tend,steps);
toc()

%% show the animation on screen
figure(1)
for t_ii = 1:length(t)
  FEMtrimesh(FEMmesh,u_dyn(:,t_ii))
  xlabel('x'); ylabel('y')
  caxis([-0.05 0.05]);
  axis([-R R -R R -0.05 0.05])
  text(4,-2,0.04,sprintf('t=%2.1f',t(t_ii)))
  drawnow();
  pause(0.3)
endfor


