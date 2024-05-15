## -*- texinfo -*-
## @deftypefn  {} {} BendingBeam1D.m
##
## This is a demo file  inside the `doc/Examples/Beam/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn
clear *
L = 3; EI = 1;
CASE = 4; %% select:    single run (1); parametrized (2)
          %% BVP1DNL(): single run (3); parametrized (4)
switch CASE
case 1
  F2 = 1.5; %% try values of 0.1 0.5 1.0 and 2
  N = 1000; s = linspace(0,L,N)';
  [sn,alpha] = BVP1D(s,1,0,0,1,F2/EI,0,[0,0]);
  figure(1); plot(sn,alpha); xlabel('arclength s'); ylabel('angle \alpha')
  xGauss = FEM1DGaussPoints(sn);
  [xGauss,Nodes2GaussU] = FEM1DGaussPoints(sn);
  for jj = 1:20
    [dalpha,ddalpha] = FEM1DEvaluateDu(sn,alpha);  %% evaluate derivative at nodes
    RHS = ddalpha + F2/EI*cos(alpha);
    alphaGauss = Nodes2GaussU*alpha;               %% evaluate u at Gauss points
    [sn,phi] = BVP1D(s,1,0,+F2/EI*sin(alphaGauss),1,RHS,0,[0,0]);
    disp(sprintf('max(abs(phi)) = %g , max(abs(RHS)) = %g',max(abs(phi)), max(abs(RHS))))
    alpha = alpha + phi;
    figure(2); plot(sn,alpha); xlabel('x'); ylabel('angle \alpha')
    pause(0.2)
  endfor
  x = cumtrapz(sn,cos(alpha)); y = cumtrapz(sn,sin(alpha));
  figure(3); plot(x,y); xlabel('x'); ylabel('y')
case 2
  F2_List = [0.25:0.25:2];
  N = 100; s = linspace(0,L,N)';
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
      %%disp(sprintf('max(abs(phi)) = %g , max(abs(RHS)) = %g',max(abs(phi)), max(abs(RHS))))
      alpha = alpha + phi;
    endfor % jj
    x = cumtrapz(sn,cos(alpha)); y = cumtrapz(sn,sin(alpha));
    figure(3); plot(x,y); xlabel('x'); ylabel('y'); hold on
  endfor %% F2
  hold off
case 3
  F2 = 1.5; %% try values of 0.1 0.25 0.5 1.0 and 2
  f = {@(s,alpha)F2/EI*cos(alpha), @(s,alpha)-F2/EI*sin(alpha)};
  N = 100; s = linspace(0,L,N)';
  %%% generate a good initial guess
  [sn,alpha0] = BVP1D(s,1,0,0,1,F2/EI,0,[0,0]);  %% as solution of linear BVP
  alpha0 = @(s)F2/(2*EI)*(L^2-(L-s).^2);           %% use the analytical solution
  [sn,alpha,inform] = BVP1DNL(s,1,0,0,1,f,0,[0,0],alpha0,
                             'tol',1e-8,'MaxIter',50,'Display','iter');
  inform
  figure(1); plot(sn,alpha); xlabel('arclength s'); ylabel('angle \alpha')
  x = cumtrapz(sn,cos(alpha)); y = cumtrapz(sn,sin(alpha));
  figure(3); plot(x,y); xlabel('x'); ylabel('y');
case 4
  F2_List = [0:0.25:2];
  N = 100; s = linspace(0,L,N)'; sn = sort([s; s(1:end-1)+diff(s)/2]);
  alpha = 0;
  for F2 = F2_List
    f = {@(s,al)F2/EI*cos(al), @(s,al)-F2/EI*sin(al)};
    [sn,alpha,inform] = BVP1DNL(s,1,0,0,1,f,0,[0,0],alpha);
    inform
  endfor
  inform
  figure(1); plot(sn,alpha*180/pi); xlabel('arclength s'); ylabel('angle \alpha [deg]')
  x = cumtrapz(sn,cos(alpha)); y = cumtrapz(sn,sin(alpha));
  figure(3); plot(x,y); xlabel('x'); ylabel('y');
endswitch

