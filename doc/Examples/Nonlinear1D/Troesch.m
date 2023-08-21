N = 51; interval = linspace(0,1,N)'; BCleft = 0; BCright = 1;
CASE = 1;
switch CASE
case 1
 alpha = 0.5;
 f = {@(x,y) -alpha*sinh(alpha*y), @(x,y)-alpha^2*cosh(alpha*y)};
 [x,y,inform] = BVP1DNL(interval,1,0,0,1,f,BCleft,BCright,@(x)x);
 figure(1); plot(x,y); xlabel('x'); ylabel('y(x)')
 xd = [ 0.1:0.1:0.9]'; yd = pwquadinterp(x,y,xd);
 Results = [xd,yd]

case 2
 figure(1); clf; hold on; box on
 y = 1;
 for alpha = [0.792, 1.151 1.753 2.394 3.308 4.129 5.0]
  f = {@(x,y) -alpha*sinh(alpha*y), @(x,y)-alpha^2*cosh(alpha*y)};
  [x,y] = BVP1DNL(interval,1,0,0,1,f,BCleft,BCright,y);
  figure(1); plot(x,y); xlabel('x'); ylabel('y(x)'); drawnow()
  [y0,dy0] = pwquadinterp(x,y,0);
  disp(sprintf("alpha = %#6.5g, y'(0)= %g",alpha,dy0))
 endfor

case 3
 y = 1;
 figure(1); clf; hold on; box on
 for alpha = 0:0.5:8.5
  f ={@(x,y) -alpha*sinh(alpha*y), @(x,y)-alpha^2*cosh(alpha*y)};
  [x,y] = BVP1DNL(interval,1,0,0,1,f,BCleft,BCright,y);
  figure(1); plot(x,y); xlabel('x'); ylabel('y(x)'); drawnow()
  [y0,dy0] = pwquadinterp(x,y,0);
  disp(sprintf("alpha = %#6.5g, y'(0)= %g",alpha,dy0))
 endfor
endswitch

