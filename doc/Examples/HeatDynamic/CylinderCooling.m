R = 2;  N = 50; alpha = linspace(0,2*pi*N/(N-1),N)';
Tend = 2; Nt = 60; %% number of shown time steps
FEMmesh = CreateMeshTriangle('circle',[R*cos(alpha),1.5*R*sin(alpha),-ones(size(alpha))],0.1);

figure(1); FEMtrimesh(FEMmesh)
           xlabel('x') ; ylabel('y'); axis equal
FEMmesh = MeshUpgrade(FEMmesh,'cubic');

function res = u_init(xy)
  x = xy(:,1); y = xy(:,2);
  res = exp(-(x-0.5).^2-2*y.^2) .* (2^2 -x.^2-y.^2);
endfunction

figure(2); FEMtrimesh(FEMmesh,u_init(FEMmesh.nodes))
           xlabel('x'); ylabel('y'); zlabel('temperature');

[u,t] = IBVP2Dsym(FEMmesh,1,1,0, 0, 0,  0,  0, 'u_init',0, Tend, [Nt,10]);

if 1  %% show figures
figure(3); FEMtrimesh(FEMmesh,u(:,Nt/4+1))
           xlabel('x'); ylabel('y'); zlim([0,1]); view([10 30]); caxis([0,1])
	   text(-1.8,-2,0.9,sprintf('t = %4.2f',t(Nt/4+1))); zlabel('temperature')
figure(4); FEMtrimesh(FEMmesh,u(:,Nt/2+1))
           xlabel('x'); ylabel('y'); zlim([0,1]); view([10 30]); caxis([0,1])
	   text(-1.8,-2,0.9,sprintf('t = %4.2f',t(Nt/2+1))); zlabel('temperature')
figure(5); FEMtrimesh(FEMmesh,u(:,3*Nt/4+1))
           xlabel('x'); ylabel('y'); zlim([0,1]); view([10 30]); caxis([0,1]);
      	   text(-1.8,-2,0.9,sprintf('t = %4.2f',t(3*Nt/4+1))); zlabel('temperature')
endif

if 0  %% an animation
  figure(11)
  steps = 2;
  for jj = 0:30
    FEMtrimesh(FEMmesh,u(:,jj*steps+1))
    text(-1.8,-2,0.9,sprintf('t = %4.2f',t(jj*steps+1))); zlabel('temperature')
    xlabel('x'); ylabel('y'); zlim([0,1]); view([10 30]); caxis([0,1])
    pause(0.2)
  endfor
endif

if 0 %% decay at center
  x = linspace(-R,R,31);
  u_center = zeros(length(x),length(t));
  for jj = 1:length(t)
    u_center(:,jj) = FEMgriddata(FEMmesh,u(:,jj),x,zeros(size(x)));
  endfor

  figure(21); mesh(t,x,u_center)
              xlabel('time t'); ylabel('x'); zlabel('temperatur')
  figure(22); contour(t,x,u_center,51)
              xlabel('time t'); ylabel('x');
  figure(23); plot(t,u_center(16,:))
              xlabel('time t'); ylabel('temperature at center')
  figure(24); semilogy(t,u_center(16,:))
              xlabel('time t'); ylabel('temperature at center')

  p = polyfit(t(40:end),log(u_center(16,40:end)),1)
  EigVal = BVP2Deig(FEMmesh,1,0,1,0,3)
endif
