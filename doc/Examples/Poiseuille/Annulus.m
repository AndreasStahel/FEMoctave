%% Poisseuille flow in an annulus
G = 1; mu = 1; R1 = 1; R2 = 2;
Interval = linspace(R1,R2,21)';
[r,u] = BVP1D(Interval,@(r)r,0,0,1,@(r)G/mu*r,0,0);

figure(1); plot(r,u); xlabel( 'radius r'); ylabel('velocity u')

u_exact = G/(4*mu)*((R1^2-r.^2)+(R2^2-R1^2)*ln(r/R1)/ln(R2/R1));
figure(2); plot(r,u-u_exact); xlabel( 'radius r'); ylabel('u-u_{exact}')

Flow = FEM1DIntegrate(r,2*pi*r.*u)

