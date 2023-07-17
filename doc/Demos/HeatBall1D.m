R = 5;  f = +3;
N = 10; r = linspace(0,R,N);
[r,u] = BVP1D(r,@(r)r.^2,0,0,@(r)r.^2,f,[0,0],0);
figure(1); plot(r,u)
           xlabel('radius r'); ylabel('temperature u')

r_fine = linspace(0,R,1001);
[u_fine,du_fine] = pwquadinterp(r,u,r_fine);

figure(2); plot(r_fine,f*4/3*pi*r_fine.^3,r_fine,-4*pi*r_fine.^2.*du_fine)
           xlabel('radius r'); ylabel('thermal energy')
           legend('energy generated','energy flux','location','northwest')

