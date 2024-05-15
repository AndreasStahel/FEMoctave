## -*- texinfo -*-
## @deftypefn  {} {} SchroedingerHarmonic.m
##
## This is a demo file  inside the `doc/Examples/Schroedinger/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

x_max = 6; interval = linspace(-x_max,x_max,100)';
BCleft = 0; BCright = 0;
[x,eVal, eVec] = BVP1Deig(interval,1,0,@(x)x.^2,1,BCleft,BCright,6);

Eigenvalues = eVal'
figure(1); plot(x,eVec(:,1:3)); xlabel('x'); ylabel('u'); xlim([-x_max,+x_max])
           legend('1','2','3')
figure(2); plot(x,eVec(:,4:6)); xlabel('x'); ylabel('u'); xlim([-x_max,+x_max])
          legend('4','5','6')

figure(11); plot(x,eVec(:,1:3).^2); xlabel('x'); ylabel('u'); xlim([-x_max,+x_max])
           legend('1','2','3'); xlim(3.5*[-1,1])
figure(12); plot(x,eVec(:,4:6).^2); xlabel('x'); ylabel('u'); xlim([-x_max,+x_max])
          legend('4','5','6'); xlim(3.5*[-1,1])
