## -*- texinfo -*-
## @deftypefn  {} {} Test_C1Conforming.m
##
## This is a demo file  inside the `doc/Demos/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

clear *
CASE = 3;

switch CASE
case 1
 n = 8; x = linspace(-1,1,n+1)';
 [x,u] = BVP1D(x,1,0,0,@(x)sign(x),1,0,0);
 figure(1); plot(x,u,'+-')
            xlabel('x'); ylabel('u')

 x_fine = linspace(-1,1,10001)'; [u_fine,du_fine] = pwquadinterp(x,u,x_fine);
 du = FEM1DEvaluateDu(x,u);
 figure(2); plot(x_fine,du_fine,x,du,'+')
            xlabel('x'); ylabel('du/dx')
            legend('interpolated','at nodes')
case 2
 n = 8; x = linspace(-1,1,n+1)';
 [x,u] = BVP1D(x,1,@(x)sign(x-0.5),0,1,@(x)exp(x),0,0);
 figure(1); plot(x,u)
            xlabel('x'); ylabel('u')

 x_fine = linspace(-1,1,10001)'; [u_fine,du_fine] = pwquadinterp(x,u,x_fine);
 du = FEM1DEvaluateDu(x,u);
 figure(2); plot(x_fine,du_fine,x,du,'+')
            xlabel('x'); ylabel('du/dx')
            legend('interpolated','at nodes')
case 3
 n = 2; x = linspace(-1,1,n+1)';
 [x,u] = BVP1D(x,1,0,0,1,@(x)exp(x),exp(-1),exp(1));
 figure(1); plot(x,u)
            xlabel('x'); ylabel('u')

 x_fine = linspace(-1,1,10001)'; [u_fine,du_fine] = pwquadinterp(x,u,x_fine);
 du = FEM1DEvaluateDu(x,u);
 figure(2); plot(x_fine,du_fine,x,du,'+')
            xlabel('x'); ylabel('du/dx')
            legend('interpolated','at nodes')

figure(3); plot(x,u,x_fine,u_fine,x_fine,exp(x_fine))
           xlabel('x'); ylabel('u')
           legend('at nodes','interpolated','exact','location','northwest')

figure(4); plot(x,u-exp(x),'+',x_fine,u_fine-exp(x_fine))
           xlabel('x'); ylabel('error of u')
           legend('at nodes','interpolated','location','northwest')

figure(5); plot(x,du-exp(x),'+',x_fine,du_fine-exp(x_fine))
           xlabel('x'); ylabel('error of du/dx')
           legend('at nodes','interpolated','location','northwest')
endswitch
