clear *
L = 3; EI = 1;
%% Newton's method
flag = 1;  %% select: single run (1); parametrized (2)
switch flag
case 1
  F2 = 0.5; %% try values of 0.1 0.5 1.5 and 2
  N = 100; s = linspace(0,L,N);
  [sn,alpha] = BVP1D(s,1,0,0,1,F2/EI,0,[0,0]);
  figure(1); plot(sn,alpha); xlabel('arclength s'); ylabel('angle \alpha')
  xGauss = FEM1DGaussPoints(sn);
  for jj = 1:20
    [dalpha,ddalpha] = FEM1DEvaluateDu(sn,alpha);
    RHS = ddalpha + F2/EI*cos(alpha);
    alphaGauss = pwquadinterp(sn,alpha,xGauss); %% evaluate u at Gauss points
    [sn,phi] = BVP1D(s,1,0,F2/EI*sin(alphaGauss),1,RHS,0,[0,0]);
    maxPHI_RHS= [max(abs(phi)) max(abs(RHS))]
    alpha = alpha + phi;
    figure(2); plot(sn,alpha); xlabel('x'); ylabel('angle \alpha')
    pause(0.5)
  endfor
  x = cumtrapz(sn,cos(alpha)); y = cumtrapz(sn,sin(alpha));
  figure(3); plot(x,y); xlabel('x'); ylabel('y')
case 2
  F2_List = [0.25:0.25:2];
  N = 100; s = linspace(0,L,N);
  [sn,alpha] = BVP1D(s,1,0,0,1,F2_List(1)/EI,0,[0,0]);
  figure(1); plot(sn,alpha); xlabel('arclength s'); ylabel('angle \alpha')
  xGauss = FEM1DGaussPoints(sn);
  figure(3); clf; hold off
  for F2 = F2_List
    for jj = 1:10
      [dalpha,ddalpha] = FEM1DEvaluateDu(sn,alpha);
      RHS = ddalpha + F2/EI*cos(alpha);
      alphaGauss = pwquadinterp(sn,alpha,xGauss); %% evaluate u at Gauss points
      [sn,phi] = BVP1D(s,1,0,F2/EI*sin(alphaGauss),1,RHS,0,[0,0]);
      alpha = alpha + phi;
    endfor % jj
    x = cumtrapz(sn,cos(alpha)); y = cumtrapz(sn,sin(alpha));
    figure(3); plot(x,y); xlabel('x'); ylabel('y'); hold on
  endfor %% F2
  hold off
endswitch

