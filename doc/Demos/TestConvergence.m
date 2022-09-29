clear *
  function res = u_exact(xy)  res =   sin(xy(:,1)).*sin(xy(:,2));   endfunction
  function res = f(xy)        res = 2*sin(xy(:,1)).*sin(xy(:,2));   endfunction
  function res = u_y(xy)      res =   sin(xy(:,1)).*cos(xy(:,2));   endfunction

a = 1; b0 = 0;  gN2 = 0;
N = 6;
Npow = 6;  % use Npow = 6 for final run
for ii = 1:Npow
  Ni = N*2^(ii-1);   h(ii) = 1/(Ni); area = 0.5/(Ni)^2;
  FEMmesh1 = CreateMeshTriangle('TestConvergence',[0 0 -1; 1 0 -1; 1 1 -2;0 1 -1],area);
  FEMmesh2 = CreateMeshTriangle('TestConvergence',[0 0 -1; 1 0 -1; 1 1 -2;0 1 -1],4*area);
  FEMmesh2 = MeshUpgrade(FEMmesh2,'quadratic');
  FEMmesh3 = CreateMeshTriangle('TestConvergence',[0 0 -1; 1 0 -1; 1 1 -2;0 1 -1],9*area);
  FEMmesh3 = MeshUpgrade(FEMmesh3,'cubic');

%%% solve with first order elements
  u1 = BVP2Dsym(FEMmesh1,a,b0,'f','u_exact','u_y',gN2);
  Difference(ii) = sqrt(FEMIntegrate(FEMmesh1,(u1-u_exact(FEMmesh1.nodes)).^2));
  [ux,uy] = FEMEvaluateGradient(FEMmesh1,u1);
  DifferenceUy(ii) = sqrt(FEMIntegrate(FEMmesh1,(uy-u_y(FEMmesh1.nodes)).^2));

%%% now for second order elements
  u2 = BVP2Dsym(FEMmesh2,a,b0,'f','u_exact','u_y',gN2);
  DifferenceQ(ii) = sqrt(FEMIntegrate(FEMmesh2,(u2-u_exact(FEMmesh2.nodes)).^2));
  [ux,uy] = FEMEvaluateGradient(FEMmesh2,u2);
  DifferenceUyQ(ii) = sqrt(FEMIntegrate(FEMmesh2,(uy-u_y(FEMmesh2.nodes)).^2));

%%% now for third order elements
  u3 = BVP2Dsym(FEMmesh3,a,b0,'f','u_exact','u_y',gN2);
  DifferenceC(ii) = sqrt(FEMIntegrate(FEMmesh3,(u3-u_exact(FEMmesh3.nodes)).^2));
  [ux,uy] = FEMEvaluateGradient(FEMmesh3,u3);
  DifferenceUyC(ii) = sqrt(FEMIntegrate(FEMmesh3,(uy-u_y(FEMmesh3.nodes)).^2));
endfor

figure(1)
plot(log10(h),log10(Difference), '+-',log10(h),log10(DifferenceUy), '+-',
     log10(h),log10(DifferenceQ),'+-',log10(h),log10(DifferenceUyQ),'+-',
     log10(h),log10(DifferenceC),'+-',log10(h),log10(DifferenceUyC),'+-')
xlabel('log_{10}(h)'); ylabel('log_{10}(difference)')
legend('linear, u-u_e','linear, d/dy (u-u_e)',
       'quad, u-u_e','quad, d/dy (u-u_e)','cubic, u-u_e','cubic, d/dy (u-u_e)',
       'location','southeast')
xlim([-2.5,-0.5])

