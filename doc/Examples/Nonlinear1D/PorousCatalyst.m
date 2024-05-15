## -*- texinfo -*-
## @deftypefn  {} {} PorousCatalyst.m
##
## This is a demo file  inside the `doc/Examples/Nonlinear1D/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

N = 51; interval = linspace(0,1,N)';
a  = 2; gamma = 20; beta = 0.05;
BCleft = [0,0]; BCright = 1;
CASE = 1; %%% 1: one run, 2: parametric, 3: second parameter set, 4: a third parameter set
switch CASE
case 1
alpha = 1;
f = {@(x,y) -alpha^2*y.*exp(gamma*beta*(1-y)./(1+beta*(1-y))),
     @(x,y) -alpha^2*exp(gamma*beta*(1-y)./(1+beta*(1-y))).*(1-y.*(gamma*beta./((1+beta*(1-y)).^2)))};
[x,y] = BVP1DNL(interval,1,@(x)-a./x,0,1,f,BCleft,BCright,1);
figure(1); plot(x,y); xlabel('x'); ylabel('y'); box on

case 2
y = 1;
figure(2); clf; hold on; box on; xlabel('x'); ylabel('y')
for alpha = 0:0.2:1
  f = {@(x,y) -alpha^2*y.*exp(gamma*beta*(1-y)./(1+beta*(1-y))),
       @(x,y) -alpha^2   *exp(gamma*beta*(1-y)./(1+beta*(1-y))).*...
                     (1-gamma*beta*y./( (1+beta*(1-y)).^2) )};
  [x,y] = BVP1DNL(interval,1,@(x)-a./x,0,1,f,BCleft,BCright,y);
  plot(x,y)
  disp(sprintf('alpha = %#3.2f, y(0) = %#g',alpha,y(1)))
endfor
legend('\alpha = 0.00','\alpha = 0.20','\alpha = 0.40','\alpha = 0.60',
       '\alpha = 0.80','\alpha = 1.00','location','southeast')
figure(1); plot(x,y); xlabel('x'); ylabel('y'); box on

case 3
a = 0; gamma = 20; beta = 0.4;
y = 1;
figure(2); clf; hold on; box on; xlabel('x'); ylabel('y')
for alpha = [0:0.05:0.35,0.36:0.001:0.370]
  f = {@(x,y) -alpha^2*y.*exp(gamma*beta*(1-y)./(1+beta*(1-y))),
       @(x,y) -alpha^2   *exp(gamma*beta*(1-y)./(1+beta*(1-y))).*(1-gamma*beta*y./( (1+beta*(1-y)).^2) )};
  [x,y] = BVP1DNL(interval,1,@(x)-a./x,0,1,f,BCleft,BCright,y);
  plot(x,y)
  disp(sprintf('alpha = %#4.3f, y(0) = %#g',alpha,y(1)))
endfor
figure(1); plot(x,y); xlabel('x'); ylabel('y'); box on

case 4
a = 2; gamma = 20; beta = 0.2;
y = 1;
figure(2); clf; hold on; box on; xlabel('x'); ylabel('y')
for alpha = [0:0.1:1.6, 1.65:0.05:2]
  f = {@(x,y) -alpha^2*y.*exp(gamma*beta*(1-y)./(1+beta*(1-y))),
       @(x,y) -alpha^2   *exp(gamma*beta*(1-y)./(1+beta*(1-y))).*(1-gamma*beta*y./( (1+beta*(1-y)).^2) )};
  [x,y] = BVP1DNL(interval,1,@(x)-a./x,0,1,f,BCleft,BCright,y);
  plot(x,y)
  %%disp(sprintf('alpha = %#4.3f, y(0) = %#g',alpha,y(1)));
endfor
figure(1); plot(x,y); xlabel('x'); ylabel('y'); box on
x_n = [0:0.2:0.8]; y_n = pwquadinterp(x,y,x_n);
disp(sprintf('y(%g)=%5.4e, y(%g)=%5.4e, y(%g)=%5.4e, y(%g)=%5.4e, y(%g)=%5.4e',
x_n(1),y_n(1),x_n(2),y_n(2),x_n(3),y_n(3),x_n(4),y_n(4),x_n(5),y_n(5)))
endswitch



