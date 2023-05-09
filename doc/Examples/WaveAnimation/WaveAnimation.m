if 0 %% linear elements
  FEMmesh = CreateMeshRect(linspace(0,pi,101),linspace(-pi,pi,101),-1,-2,-2,-2);
else %% quadratic elements
  FEMmesh = CreateMeshRect(linspace(0,pi,51),linspace(-pi,pi,51),-1,-2,-2,-2);
  FEMmesh = MeshUpgrade(FEMmesh);
endif
x = FEMmesh.nodes(:,1); y = FEMmesh.nodes(:,2);

m=1; alpha=0.0; a=1; b0=0; bx=0; by=0; f=0; 
gD=0; gN1=0; gN2=0;
t0=0; tend=3 ; steps = [150,10];

u0 = exp(-25*((x-1).^2+(y-0).^2)); 
v0 = zeros(length(FEMmesh.nodes),1);
tic()
[u_dyn,t] = I2BVP2D(FEMmesh,m,alpha,a,b0,bx,by,f,gD,gN1,gN2,u0,v0,t0,tend,steps);
toc()

figure(1)   % show animation
for t_ii = 1:length(t)
  FEMtrimesh(FEMmesh,u_dyn(:,t_ii))
  axis([0 pi -pi pi -0.2 0.4])
  xlabel('x'); ylabel('y')
  drawnow(); %%pause(0.1)
endfor
