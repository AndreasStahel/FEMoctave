## -*- texinfo -*-
## @deftypefn  {} {} SoundWaveSpherical.m
##
## This is a demo file  inside the `doc/Examples/SoundWave/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

R = 2; H = 2; N = 60;
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
[u,t] = I2BVP2D(FEMmesh,'rho',0,'rho',0,0,0,0,0,0,0,'u_0','v_0',0,1.75,[100,10]);
ComputationTime = toc()

figure(1); clf
if 0 %% animation
  for jj = 1:length(t)
    FEMtrimesh(FEMmesh,u(:,jj))
    xlabel('rho'); ylabel('z');
    zlim([-0.5 0.5]); caxis(0.1*[-0.5,0.5])
    pause(0.1)
  endfor
else
  FEMtrimesh(FEMmesh,u(:,end))
  xlabel('\rho'); ylabel('z')
endif

max_u = max(u) - min(u); t_start = find(t>0.6,1);
t_tail = t(t_start:end)';

[p,~,~,p_var] = LinearRegression(1./t_tail,max_u(t_start:end)');
figure(11); plot(t,max_u,t_tail, p./t_tail)
            xlabel('time t'); ylabel ('amplitude'); xlim([0,max(t)])
            legend('amplitude','best fit')

