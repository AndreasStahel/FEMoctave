rho = 1100; E = 1e6; nu = 0.47;  %% Silicone rubber
H = 0.1; R = 0.2; W = 0.01;
Contour = [0 H -11;0 H-2*W -22; R-2*W H-2*W -22; R-W H-2*W -22;R-W 0 -21;
           R 0 -22; R H-W -22; R-W H -22];
Mesh = CreateMeshTriangle('RubberBox',Contour,3e-5);
Mesh = MeshUpgrade(Mesh,'cubic');

function res=fr(xy,dummy)
  freq = 10; omega = freq*2*pi; rho = 1100;
  res = rho*xy(:,1)*omega^2;
endfunction

[ur,uz] = AxiStress(Mesh,E,nu,{'fr',0},{0,0},{0,0});

factor = 1;
figure(10); trimesh(Mesh.elem,Mesh.nodes(:,1)+factor*ur,...
                   Mesh.nodes(:,2)+factor*uz,'color','red','linewidth',2);
hold on ;  trimesh(Mesh.elem,Mesh.nodes(:,1),Mesh.nodes(:,2),'color','green','linewidth',1);
           hold off; axis equal; xlabel('x'); ylabel('y');
figure(11); FEMtrimesh(Mesh,ur)
            xlabel('r'); ylabel('z'); zlabel('u_r')
figure(12); FEMtrimesh(Mesh,uz)
            xlabel('r'); ylabel('z'); zlabel('u_z')

[sigma_x,sigma_y,sigma_z,tau_xz] = EvaluateStressAxi(Mesh,ur,uz,E,nu);
vonMises = EvaluateVonMisesAxi(sigma_x,sigma_y,sigma_z,tau_xz);

figure(13); FEMtrimesh(Mesh,vonMises/1e6)
            xlabel('r'); ylabel('z'); zlabel('von Mises [MPa]'); view([35 30])
figure(14); clf; FEMtricontour(Mesh,vonMises/1e6)
            xlabel('r'); ylabel('z'); zlabel('von Mises [MPa]')
            hold on; plot([Contour(:,1);Contour(1,1)],[Contour(:,2);Contour(1,2)],'k')
            hold off; axis equal; colorbar; %%shading interp
            title('von Mises stress [MPa]')
