%% solve the dynamic Schroedinger equation for a bouncing ball
Lm = -10; Lp = +10; BCleft = 0; BCright = 0; w = 1; a = 1; b = 0; d = 1; f = 0;
V0 = 21; Width = 1 ; speed = -4;
c = @(x) 15*x;
interval = linspace(Lm,Lp,1001)'; t0 = 0; t_end = 3.5;  steps = [400,50];
sigma = sqrt(2); Position = -5 ;
u0 = @(x)exp(-0.5*sigma^(-2)*(x+Position).^2).^2; cu = 1/sqrt(quad(u0,Lm,Lp));
u0 = @(x)cu*exp(-0.5*sigma^(-2)*(x+Position).^2).*exp(i*speed/2*x); %% nonzero initial velocity


if 1  %% dynamic, complex
  solver = 'RK';
  [x,u,t] = IBVP1D(interval,-w*i,a,b,c,d,f,BCleft,BCright,u0,t0,t_end,steps,'solver',solver);
  figure(11); mesh(t,x,real(u)); xlabel('time t'); ylabel('position x'); xlim([0,t_end])
              title('real(u)')
  figure(12); mesh(t,x,imag(u)); xlabel('time t'); ylabel('position x'); xlim([0,t_end])
              title('imag(u)')

  figure(13); mesh(t,x,abs(u).^2); xlabel('time t'); ylabel('position x'); xlim([0,t_end])
              title('|u|^2')

  figure(14); contour(t,x,abs(u).^2); xlabel('time t'); ylabel('position x')
              title('contours of |u|^2')
  figure(15); plot(x,real(u(:,end)),'b',x,imag(u(:,end)),'g',x,abs(u(:,end)).^2,'r')
              xlabel('position x'); title(['at t = ',num2str(t_end)])
              legend('real(u)','imag(u)','|u|^2')
endif

%% u_norms = sum(abs(u).^2)*(x(2)-x(1)); %% is OK, i.e. all = 1

printing = 0;
if printing
  figure(11); print -dpng SchroedingerBouncingBallReal.png
  figure(12); print -dpng SchroedingerBouncingBallImag.png
  figure(13); print -dpng SchroedingerBouncingBallNorm.png
  figure(14); print -dpdfcrop SchroedingerBouncingBallContour.pdf
  figure(15); print -dpdfcrop SchroedingerBouncingBallFinal.pdf
endif

