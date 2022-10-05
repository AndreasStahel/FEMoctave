l = 1; h = 0.1; L = 4; d = 2; H = 2;
FEMmesh = CreateMeshTriangle('test',[-l,0,-2;0 0 -2;0,-d,-2;L,-d,-2; L,H,-2;0,H,-2;0,h,-2;-l,h,-2],0.01);
figure(1); FEMtrimesh(FEMmesh)
           xlabel('x'); ylabel('y'); axis equal
FEMmesh = MeshUpgrade(FEMmesh,'cubic');

function res = f(xy,t)
  res = cos(0.2*pi*t)*ones(size(xy,1),1);
  res(xy(:,1)>-0.9) = 0;
endfunction

function res = u0(xy)   res = zeros(size(xy,1),1); endfunction

m = 1; a = 1; b0 = 0; bx = by = 0; f = 0; gn1 = gn2 = 0;
tic();
%%[u,t] = IBVP2D(FEMmesh,m,a,b0,bx,by,'f',0,gn1,gn2,'u0',0,30,[2*60,10]);
[u,t] = IBVP2Dsym(FEMmesh,m,a,b0,'f',0,gn1,gn2,'u0',0,30,[2*60,10]);
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
    xlim([0,L])
%    view(0,90); xlim([0,L]);  ylim([-d,H]);
    pause(0.1);
  endfor
else
  FEMtrimesh(FEMmesh,u(:,end))
  xlabel('x'); ylabel('y')
  zlim(umax*[-1 1]); caxis(0.3*umax*[-1 1]);
  text(0.8*L,0.8*H,umax,sprintf('t = %4.2f',t(end)))
%%  view(0,90); xlim([0,L]);  ylim([-d,H]);
endif

x = linspace(0,L,101);
u_line = zeros(size(t,1),size(x,2));
for jj = 1:length(t)
  u_line(jj,:) = FEMgriddata(FEMmesh,u(:,jj),x,ones(size(x)));
endfor

figure(4); mesh(x,t,u_line)
           xlabel('x'); ylabel('t');
printing = 0;
if printing
  print -dpng ../../doc/HeatWaveSlice.png
endif

figure(5); contour(x,t,u_line,0.003*[-1:0.1:+1])
           xlabel('x'); ylabel('t');
if printing
  print -dpdfcrop ../../doc/HeatWaveContour.pdf
endif


