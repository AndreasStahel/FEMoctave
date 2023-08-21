CASE = 1;  %% use 1 for the codes successive substitution, 2 for BVP1DNL()
switch CASE
  case 1
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

  figure(1); plot(phi/pi*180,90-theta/pi*180);
             xlabel('\phi [deg]'); ylabel('90-\theta [deg]'); hold on
  x = sin(theta).*cos(phi); y = sin(theta).*sin(phi); z = cos(theta);
  Dtheta = FEM1DEvaluateDu(phi,theta);
  figure(2); plot3(x,y,z);
             xlabel('x'); ylabel('y'); zlabel('z'); view([-20,10]); hold on
  disp(sprintf("Iteration 1: L = %f km",R*trapz(phi,sqrt(sin(theta).^2+Dtheta.^2))))
  figure(3); plot(phi*180/pi,[x,y,z]*n);
             xlabel('\phi [deg]'); hold on;
  phiG = FEM1DGaussPoints(phi);
  for jj = 2:5;
    [thetaG,DthetaG] = pwquadinterp(phi,theta,phiG);
    a = (sin(thetaG).^2 + DthetaG.^2).^(-1/2);
    d = -sin(thetaG).*cos(thetaG)./sqrt(sin(thetaG).^2 + DthetaG.^2);
    [phi,theta] = BVP1D(interval,a,0,0,d,1,BCleft,BCright);
    Dtheta = FEM1DEvaluateDu(phi,theta);
    disp(sprintf("Iteration %i: L = %f km",
                jj,R*trapz(phi,sqrt(sin(theta).^2+Dtheta.^2))))
    figure(1); plot(phi/pi*180,90-theta/pi*180);
    x = sin(theta).*cos(phi); y = sin(theta).*sin(phi); z = cos(theta);
    figure(2); plot3(x,y,z);
    figure(3); plot(phi*180/pi,[x,y,z]*n);
    pause(0.2)
  endfor %% jj
  figure(1); legend('1','2','3','4','5'); hold off;
  figure(2); hold off; figure(3); hold off
case 2
  N = 51; R = 6300;
  Angles_ZH = [8,90-47]/180*pi;   Angles_SLC = [-112,90-41]/180*pi;
  BCleft = Angles_SLC(2); BCright = Angles_ZH(2);
  interval = linspace(Angles_SLC(1),Angles_ZH(1),N);
  a =  @(x,u,du)(sin(u).^2+du.^2).^(-0.5);
  f = {@(x,u,du)(-cos(u).*sin(u))./sqrt(sin(u).^2+du.^2),
       @(x,u,du)-(1-2*sin(u).^2)./sqrt(sin(u).^2+du.^2)...
                +(cos(u).^2.*sin(u).^2)./(sqrt(sin(u).^2+du.^2).^3),
       @(x,u,du)(cos(u).*sin(u).*du)./(sqrt(sin(u).^2+du.^2).^3)};
  [phi2,theta2] = BVP1D(interval,1,0,0,1,0,BCleft,BCright);  %% generate an intial guess
  [phi2,theta2,inform] = BVP1DNL(interval,a,0,0,1,f,BCleft,BCright,theta2,
                                'MaxIter',20,'Display','iter');
  figure(4); plot(phi2/pi*180,90-theta2/pi*180);
             xlabel('\phi [deg]'); ylabel('90-\theta [deg]');
  Dtheta2 = FEM1DEvaluateDu(phi2,theta2);
  disp(sprintf("BVP1DNL(): L = %f km",R*trapz(phi2,sqrt(sin(theta2).^2+Dtheta2.^2))))
endswitch


