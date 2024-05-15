## -*- texinfo -*-
## @deftypefn  {} {} WasherSpring.m
##
## This is a demo file  inside the `doc/Examples/Elasticity/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

rho = 2700; E = 70e9; nu = 0.33;  %% Aluminum
H = 0.001; Ri = 0.002; Ro = 0.005; D = 0.0006; H = 0.0004;
global Offset
Offset = 1*1e-4;

if 1 %% free sides
  Contour = [Ri H+D -22; Ri H -22;Ro-D 0 -21; Ro 0 -22; Ro D -22;Ri+D H+D -21];
elseif 1 %% clamped on the outside
  Contour = [Ri H+D -22; Ri H -22;Ro-D 0 -21; Ro 0 -12; Ro D -22;Ri+D H+D -21];
else  %% clamped on both sides
  Contour = [Ri H+D -12; Ri H -22;Ro-D 0 -21; Ro 0 -12; Ro D -22;Ri+D H+D -21];
endif

Mesh = CreateMeshTriangle('Washer',Contour,2.5e-9);
%%Mesh = MeshUpgrade(Mesh,'quadratic');
Mesh = MeshUpgrade(Mesh,'cubic');

function res = gDz(xy,dummy)
  global Offset
  res = -Offset*(xy(:,2)>Offset);
endfunction

[ur,uz] = AxiStress(Mesh,E,nu,{0,0},{0,'gDz'},{0,0});

figure(10); ShowDeformation(Mesh,ur,uz,1);
            axis equal; xlabel('r'); ylabel('z'); xticks([2:5]/1000); yticks([0:0.5:1]/1000)
figure(11); FEMtrimesh(Mesh,ur)
            xlabel('r'); ylabel('z'); zlabel('u_r')
            xticks([2:5]/1000); yticks([0:0.5:1]/1000)
Cx = [Contour(:,1);Contour(1,1)]; Cy = [Contour(:,2);Contour(1,2)];
figure(21); clf; FEMtricontour(Mesh,ur)
            hold on ;  plot(Cx,Cy,'k'); hold off
            xlabel('r'); ylabel('z'); title('u_r');
            axis equal; colorbar;xticks([2:5]/1000); yticks([0:0.5:1]/1000)

figure(12); FEMtrimesh(Mesh,uz)
            xlabel('r'); ylabel('z'); zlabel('u_z');xticks([2:5]/1000); yticks([0:0.5:1]/1000)
figure(22); clf; FEMtricontour(Mesh,uz)
            hold on ;  plot(Cx,Cy,'k'); hold off
            xlabel('r'); ylabel('z'); title('u_z');
            axis equal; colorbar;xticks([2:5]/1000); yticks([0:0.5:1]/1000)

[sigma_x,sigma_y,sigma_z,tau_xz] = EvaluateStressAxi(Mesh,ur,uz,E,nu);
figure(13); FEMtrimesh(Mesh,sigma_z*1e-6)
            xlabel('r'); ylabel('z'); zlabel('\sigma_z [MPa]');xticks([2:5]/1000); yticks([0:0.5:1]/1000)
figure(23); clf; FEMtricontour(Mesh,sigma_z/1e6,[-20:1:20]*100)
            hold on ;  plot(Cx,Cy,'k'); hold off
            xlabel('r'); ylabel('z'); title('\sigma_z [MPa]');
            axis equal; colorbar;xticks([2:5]/1000); yticks([0:0.5:1]/1000)

r = linspace(0,D,1000)';
sigma_up = FEMgriddata(Mesh,sigma_z,Ri+r,(H+D)*ones(size(r)));
figure(31); plot(Ri+r,sigma_up); xlabel('r'); ylabel('\sigma_z'); title('upper edge')
xlim([Ri,Ri+D]); xticks([2:0.2:2.6]/1000);
Force_up = 2*pi*trapz(Ri+r,sigma_up.*(Ri+r))

sigma_low = FEMgriddata(Mesh,sigma_z,Ro-D+r,zeros(size(r)));
figure(32); plot(Ro-D+r,sigma_low); xlabel('r'); ylabel('\sigma_z'); title('lower egde')
xlim([Ro-D, Ro]);  xticks([4.4:0.2:5]/1000);
Force_low = 2*pi*trapz(Ro-D+r,sigma_low.*(Ro-D+r))

s = 0.5;  %% select the height
r_mid = linspace(Ri,Ro,1000)';
sigma_mid = FEMgriddata(Mesh,sigma_z,r_mid,s*(H+D)*ones(size(r_mid)));
ind = find(isfinite(sigma_mid));
r_mid = r_mid(ind); sigma_mid = sigma_mid(ind);
figure(33); plot(r_mid,sigma_mid); xlabel('r'); ylabel('\sigma_z'); title('at half height');
            xticks([2:5]/1000);
Force_mid = 2*pi*trapz(r_mid,sigma_mid.*r_mid)



[urGP,ur_rz] = FEMEvaluateGP(Mesh,ur);
[uzGP,uz_rz] = FEMEvaluateGP(Mesh,uz);
rGP = Mesh.GP(:,1);
w = E/(2*(1+nu)*(1-2*nu))*rGP.*((1-nu)*(ur_rz(:,1).^2+uz_rz(:,2).^2+(urGP./rGP).^2)...
                                + 2*nu*(ur_rz(:,1).*uz_rz(:,2) +1./rGP.*urGP.*(ur_rz(:,1)+uz_rz(:,2))))...
      +E/(4*(1+nu))*rGP.*(ur_rz(:,2)+uz_rz(:,1)).^2;
U_elast = 2*pi*FEMIntegrate(Mesh,w);
Force_energy = 2*U_elast/Offset

