## -*- texinfo -*-
## @deftypefn  {} {} NL_Newton.m
##
## This is a demo file  inside the `doc/Demos/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

interval = linspace(-1,2,40)';
al = +0.1;
f = {@(x,u)0.5+al*(x.*u+sin(u)), @(x,u)+al*(x+cos(u))};
BCleft = 1; BCright = 4;
u0 = @(x)2+x;
[x,u,inform] = BVP1DNL(interval,1,0,0,1,f,BCleft,BCright,u0,'Display','iter');
inform
figure(1); plot(x,u); xlabel('x'); ylabel('u')
