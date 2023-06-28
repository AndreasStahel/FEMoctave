R1 = 1; R2 = 3; BCleft = 0; BCright = 1;
interval = linspace(R1,R2,11);
[r,u] = BVP1D(interval,@(r)r.^2,0,0,1,0,BCleft,BCright);
figure(4); plot(r,u,r,R2/(R2-R1)*(r-R1)./r)
           xlabel('radius r'); ylabel( 'temperature u')
           legend('u_{FEM}','u_{exact}','location','northwest')

figure(5); plot(r,u-R2/(R2-R1)*(r-R1)./r,'+-')
           xlabel('radius r'); ylabel( 'temperature u')
           legend('difference','location','southeast')

r_fine = linspace(R1,R2,501)'; u_fine = pwquadinterp(r,u,r_fine);
u_exact = R2/(R2-R1)*(r_fine-R1)./r_fine;
figure(6); plot(r_fine,u_fine-u_exact,'k',r,u-R2/(R2-R1)*(r-R1)./r,'b+')
           xlabel('radius r'); ylabel('u_{FEM}-u_{exact}')
           legend('interpolated','at nodes')
