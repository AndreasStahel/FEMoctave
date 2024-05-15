## -*- texinfo -*-
## @deftypefn  {} {} SoundWaveCylindrical.m
##
## This is a demo file  inside the `doc/Examples/SoundWave/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

R = 4.5; H = 4.5;
N = 60;
FEMmesh = CreateMeshRect(linspace(0,R,N),linspace(0,H,N),-2,-2,-2,-2);
FEMmesh = MeshUpgrade(FEMmesh,'quadratic');

function res = rho(xy,dummy); res = xy(:,1); endfunction;
function res = u_0(xy)
  r = sqrt(xy(:,1).^2+xy(:,2).^2);
  res = 1+cos(10*r);
  res(r>pi/10) = 0;
endfunction
function res = v_0(xy) ; res = zeros(size(xy,1),1); endfunction
tic();
[u,t] = I2BVP2D(FEMmesh,'rho',0,'rho',0,0,0,0,0,0,0,'u_0','v_0',0,4,[101,50]);
ComputationTime = toc()

figure(3); clf
if 0 %% animation
  for jj = 1:length(t)
    FEMtrimesh(FEMmesh,u(:,jj))
    xlabel('rho'); ylabel('z');
    zlim(0.1*[-2 2]); caxis(0.5*[-2 2])
    pause(0.1)
  endfor
else
  FEMtrimesh(FEMmesh,u(:,end))
  xlabel('\rho'); ylabel('z')
endif

if 0 %% slice along z=0
  x = linspace(0,R);
  u_slice = zeros(size(t,1),length(x));
  for jj = 1:length(t)
    u_slice(jj,:) = FEMgriddata(FEMmesh,u(:,jj),x,zeros(size(x)));
  endfor
  figure(4); mesh(x,t,u_slice)
             xlabel('position x'); ylabel('time t')
             zlim([-0.3 0.3]); caxis([-0.3 0.3]); view(75,20);
endif

max_u = max(u);
p = LinearRegression([ones(length(t(20:end)),1),1./t(20:end)'],max_u(20:end)');
figure(11); plot(t,max(u),t(20:end), p(1)+ p(2)./t(20:end))
            xlabel('time t'); ylabel('maximal value')

