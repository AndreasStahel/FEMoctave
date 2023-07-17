clear *
N = 10;              % number of elements, then 2*N+1 nodes
x = linspace(0,pi/2,N+1);
[xn,u] = BVP1D(x,1,0,0,1,@(x)-sin(x),0,1);

x_fine = linspace(0,pi/2,1001);
[u_fine du_fine, ddu_fine] = pwquadinterp(xn,u,x_fine);

figure(1); plot(x_fine,u_fine)
            xlabel('x'); ylabel('u')

figure(11); plot(x_fine,u_fine-sin(x_fine),xn,u-sin(xn),'+')
            xlabel('x'); ylabel('u_{FEM}-u_{exact}'); xlim([min(x_fine),max(x_fine)])
            legend('interpolated','at nodes')

figure(2); plot(x_fine,du_fine,x_fine,cos(x_fine))
            xlabel('x'); ylabel('du/dx'); xlim([min(x_fine),max(x_fine)])
            legend('FEM','exact')


figure(12); plot(x_fine,du_fine-cos(x_fine))
            xlabel('x'); ylabel('difference du/dx'); xlim([min(x_fine),max(x_fine)])

figure(3); plot(x_fine,ddu_fine,x_fine,-sin(x_fine))
            xlabel('x'); ylabel('d^2u/dx^2'); xlim([min(x_fine),max(x_fine)])
            legend('FEM','exact')

figure(13); plot(x_fine,ddu_fine+sin(x_fine))
            xlabel('x'); ylabel('difference d^2u/dx^2'); xlim([min(x_fine),max(x_fine)])
