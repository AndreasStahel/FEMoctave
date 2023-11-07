%% this code requires that RingVibration was run first
Mode = 7
if     Mode<16   n = floor(Mode/2)
elseif Mode<19   n = floor((Mode-1)/2)
else             n = floor((Mode - 17)/2) endif
phi = linspace(0,2*pi,200)'; x = R*cos(phi); y = R*sin(phi);
u1_r = FEMgriddata(FEMmesh,u1_all(:,Mode),x,y);
u2_r = FEMgriddata(FEMmesh,u2_all(:,Mode),x,y);
ur = cos(phi).*u1_r + sin(phi).*u2_r;
uphi = -sin(phi).*u1_r + cos(phi).*u2_r;

M = [ones(size(phi)),cos(n*phi),sin(n*phi)];  p = LinearRegression(M,ur);
dr = p(1); Amp_r = sqrt(p(2)^2+p(3)^2); phase_r = atan2(p(2),p(3))/pi*180;
disp(sprintf('Results for u_r:   dr   = %g, amplitude = %g, phase_r   = %g',dr,Amp_r,phase_r))
ur_fit = M*p;
figure(101); plot(phi*180/pi,ur,'g',phi/pi*180,ur_fit,'b'); xlim([0,360])
             xlabel('phi [deg]'); ylabel('u_r'); legend('data','fit');

p = LinearRegression(M,uphi);
dphi = p(1); Amp_phi = sqrt(p(2)^2+p(3)^2);  phase_phi = atan2(p(2),p(3))/pi*180;
disp(sprintf('Results for u_phi: dphi = %g, amplitude = %g, phase_phi = %g',dphi,Amp_phi,phase_phi))
disp(sprintf('Ratio of amplitudes = %g, phaseDiff = %g',Amp_phi/Amp_r,phase_r-phase_phi))
uphi_fit = M*p;
figure(102); plot(phi*180/pi,uphi,'g',phi/pi*180,uphi_fit,'b');  xlim([0,360])
             xlabel('phi [deg]'); ylabel('u_\phi'); legend('data','fit'); xlim([0,360])


