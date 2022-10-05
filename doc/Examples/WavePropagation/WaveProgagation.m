l = 1; h = 0.1; L = 10; d = 2; H = 9;
FEMmesh = CreateMeshTriangle('test',[-l,0,-2;0 0 -2;0,-d,-2;L,-d,-2; L,H,-2;0,H,-2;0,h,-2;-l,h,-2],0.1);
figure(1); FEMtrimesh(FEMmesh)
           xlabel('x'); ylabel('y'); axis equal
FEMmesh = MeshUpgrade(FEMmesh,'cubic');

function res = f(xy,t)
  res = sin(3*pi*t)*ones(size(xy,1),1);
  res(xy(:,1)>-0.9) = 0;
endfunction

function res = v0(xy);  res = zeros(size(xy,1),1); endfunction
function res = u0(xy)   res = zeros(size(xy,1),1); endfunction

m = 1; a = 1; b0 = 0; bx = by = 0; f = 0; gn1 = gn2 = 0;
tic();
[u,t] = I2BVP2D(FEMmesh,m,0,a,b0,bx,by,'f',0,gn1,gn2,'u0','v0',0,11,[56,10]);
SolverTime = toc()

figure(2); FEMtrimesh(FEMmesh,u(:,end))
           xlabel('x'); ylabel('y'); xlim([0,L]);

umax = 0.3*max([-min(u(:)),max(u(:))]);
figure(3)
if 0  %% animation
  for jj = 1:length(t)
    FEMtrimesh(FEMmesh,u(:,jj))
    xlabel('x'); ylabel('y')
    zlim(umax*[-1 1]); caxis(0.3*umax*[-1 1]);
    text(0.8*L,0.8*H,umax,sprintf('t = %4.2f',t(jj)))
    view(0,90); xlim([0,L]);  ylim([-d,H]);
    pause(0.1);
  endfor
else
  FEMtrimesh(FEMmesh,u(:,end))
  xlabel('x'); ylabel('y')
  zlim(umax*[-1 1]); caxis(0.3*umax*[-1 1]);
  text(0.8*L,0.8*H,umax,sprintf('t = %4.2f',t(end)))
  view(0,90); xlim([0,L]);  ylim([-d,H]);
endif

printing = 0;
if printing
  print -dpng ../../doc/WavePropagation.png
endif

