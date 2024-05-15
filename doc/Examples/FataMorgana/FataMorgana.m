## -*- texinfo -*-
## @deftypefn  {} {} FataMorgana.m
##
## This is a demo file inside the `doc/Examples/FataMorgana/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

Interval = linspace(0,1,21)';
alpha = 0.2;
f = {@(x,u,v)alpha*sqrt(1+v.^2),
     @(x,u,v)0*u;
     @(x,u,v)alpha*v./sqrt(1+v.^2)};

a =  @(x,u,v)(1+alpha*u)./sqrt(1+v.^2);

[x,u,inform] = BVP1DNL(Interval,a,0,0,1,f,0,0,0);

figure(1); plot(x,u,'b');

du = FEM1DEvaluateDu(x,u); du1 = du(1);
[x,u2] = BVP1DNL(Interval,a,0,0,1,f,0,0.05,0);
du = FEM1DEvaluateDu(x,u2); du2 = du(1);
figure(2); plot(x,u,'b',[0,1],[0,du1],'g',x,u2,'b',[0,1],[0,du2],'g')
           legend('true','visual','location','northwest')

