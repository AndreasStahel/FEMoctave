N = 51; R = 6300;
Angles_ZH = [8,90-47]/180*pi;   Angles_SLC = [-112,90-41]/180*pi;
Vec_ZH = R*[sin(Angles_ZH(2))*cos(Angles_ZH(1));
            sin(Angles_ZH(2))*sin(Angles_ZH(1));cos(Angles_ZH(2))];
Vec_SLC = R*[sin(Angles_SLC(2))*cos(Angles_SLC(1));
             sin(Angles_SLC(2))*sin(Angles_SLC(1));cos(Angles_SLC(2))];
n = cross(Vec_SLC,Vec_ZH); n = n/norm(n);
BCleft = Angles_SLC(2); BCright = Angles_ZH(2);
interval = linspace(Angles_SLC(1),Angles_ZH(1),N);
[phi,theta] = BVP1D(interval,1,0,0,1,0,BCleft,BCright);

figure(1); plot(phi/pi*180,90-theta/pi*180); xlabel('\phi [deg]'); ylabel('90-\theta [deg]'); hold on
x = sin(theta).*cos(phi); y = sin(theta).*sin(phi); z = cos(theta);
figure(2); plot3(x,y,z);
           xlabel('x'); ylabel('y'); zlabel('z'); view([-20,10]); hold on
disp(sprintf("Iteration 1: L = %f km",
              R*trapz(phi,sqrt(sin(theta).^2+Dtheta.^2))))
Dtheta = FEM1DEvaluateDu(phi,theta);
figure(3); plot(phi*180/pi,[x,y,z]*n);
           xlabel('\phi [deg]'); hold on;
phiG = FEM1DGaussPoints(phi);
for jj = 2:6;
  Dtheta = FEM1DEvaluateDu(phi,theta);
  disp(sprintf("Iteration %i: L = %f km",
              jj,R*trapz(phi,sqrt(sin(theta).^2+Dtheta.^2))))
  [thetaG,DthetaG] = pwquadinterp(phi,theta,phiG);
  a = (sin(thetaG).^2 + DthetaG.^2).^(-1/2);
  d = -sin(thetaG).*cos(thetaG)./sqrt(sin(thetaG).^2 + DthetaG.^2);
  [phi,theta] = BVP1D(interval,a,0,0,d,1,BCleft,BCright);
  figure(1); plot(phi/pi*180,90-theta/pi*180);
  x = sin(theta).*cos(phi); y = sin(theta).*sin(phi); z = cos(theta);
  figure(2); plot3(x,y,z);
  figure(3); plot(phi*180/pi,[x,y,z]*n);
  pause(1)
endfor %% jj
figure(1); legend('1','2','3','4','5','6'); hold off;
figure(2); hold off; figure(3); hold off


