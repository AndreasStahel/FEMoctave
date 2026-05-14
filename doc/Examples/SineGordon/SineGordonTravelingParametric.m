c = 0.2; L = 4;  Interval = linspace(0,L,1001);
f  = {@(x,u)sin(u),@(x,u)cos(u)};
u0 = @(x)pi+pi*sin(pi/2*x/L);
[x,u,inform] = BVP1DNL(Interval,(c^2-1),0,0,1,f,pi,[0,0],u0);
figure(1); clf; plot(x,u0(x),x,u)
           xlabel('u'); ylabel('u');
           legend('u_0','data0','location','southeast')
           hold on
c_list = [0.3 0.4 0.5 0.6 0.65 0.7 0.75 0.775 0.8];
for ii = 1:length(c_list)
  c = c_list(ii)
  [x,u,inform] = BVP1DNL(Interval,(c^2-1),0,0,1,f,pi,[0,0],u,'tol',1e-10, 'display','off');
  figure(1); plot(x,u); xlabel('u'); ylabel('u');
  pause(1)
endfor
u0 = @(x)4*atan(exp(x/sqrt(1-c^2)));
figure(2); plot(x,u,x,u0(x)); xlabel('u'); ylabel('u');
           legend('u','u_{exact}','location','northwest'); hold off

