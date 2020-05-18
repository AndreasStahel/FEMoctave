clear *
function res = u_exact(xy)  res = sin(xy(:,1)).*sin(xy(:,2));    endfunction
function res = f(xy)        res = 2*sin(xy(:,1)).*sin(xy(:,2));  endfunction
function res = u_y(xy)      res = sin(xy(:,1)).*cos(xy(:,2));    endfunction

a = 1; b0 = 0; gN2 = 0;
N = 6; %% N=6 for final run
Npow = 6;  % use Npow = 6 for final run
for ii = 1:Npow
  Ni = N*2^(ii-1)+1;   h(ii) = 1/(Ni-1); area = 0.5/(Ni-1)^2;
  tic();
%  FEMmesh = CreateMeshRect(linspace(0,1,Ni/2),linspace(0,1,Ni/2),-1,-2,-1,-1);
  FEMmesh = CreateMeshTriangle('TestConvergence',[0 0 -1; 1 0 -1; 1 1 -2;0 1 -1],4*area);
  FEMmesh = MeshUpgrade(FEMmesh);
  FEMmeshLin = MeshQuad2Linear(FEMmesh);
  disp(sprintf("\nDOF = %i and %i",FEMmesh.nDOF,FEMmeshLin.nDOF))
  time_mesh = toc()
  %%% solve with first order elements
  tic()
  u = BVP2Dsym(FEMmeshLin,a,b0,'f','u_exact','u_y',gN2);
  time_solve_lin = toc()
  Difference(ii) = sqrt(FEMIntegrate(FEMmeshLin,(u-u_exact(FEMmeshLin.nodes)).^2));
  tic()
  [ux,uy] = FEMEvaluateGradient(FEMmeshLin,u);
  time_grad_lin = toc()
  DifferenceUy(ii) = sqrt(FEMIntegrate(FEMmesh,(uy-u_y(FEMmesh.nodes)).^2));
  %%% now for second order elements
  tic();
  u = BVP2Dsym(FEMmesh,a,b0,'f','u_exact','u_y',gN2);
  time_solve_quad = toc()
  DifferenceQ(ii) = sqrt(FEMIntegrate(FEMmesh,(u-u_exact(FEMmesh.nodes)).^2));
  tic()
  [ux,uy] = FEMEvaluateGradient(FEMmesh,u);
  time_grad_quad = toc()
  DifferenceUyQ(ii) = sqrt(FEMIntegrate(FEMmesh,(uy-u_y(FEMmesh.nodes)).^2));
endfor
figure(1)
plot(log10(h),log10(Difference),'+-',log10(h),log10(DifferenceUy),'+-',
     log10(h),log10(DifferenceQ),'+-',log10(h),log10(DifferenceUyQ),'+-')
xlabel('log_{10}(h)'); ylabel('log_{10}(difference)')
legend('linear, u-u_e','linear, d/dy (u-u_e)',
       'quad, u-u_e','quad, d/dy (u-u_e)','location','southeast')
xlim([-2.5,-0.5])
