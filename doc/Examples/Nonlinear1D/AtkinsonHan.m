%% AtkiHan09 Exercise 5.4.1 p. 241
N = 501;
interval = linspace(0,1,N)';
BCleft = 1; BCright = [2,0];

u = 0;
figure(2); clf; hold on; box on
for alpha = linspace(0,1,5)
  alpha
  f = {@(x,u,du)-u.*du - alpha*u.^3 + exp(x),@(x,u,du)-du -alpha*3*u.^2 ,@(x,u,du)-u};
  [x,u,inform] = BVP1DNL(interval,1,0,0,1,f,BCleft,BCright,u,
                 'Display','iter');
  plot(x,u) ;xlabel('x'); ylabel('u(x)')
  pause(0.1)
endfor
legend('\alpha = 0.00','\alpha = 0.25','\alpha = 0.50','\alpha = 0.75','\alpha = 1.00','location','northwest')
figure(1); plot(x,u)  ; xlabel('x'); ylabel('u(x)')
