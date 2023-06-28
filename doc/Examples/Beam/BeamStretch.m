L = 3; F = 0.2; N = 10;
EA = @(x) (2-sin(x/L*pi))/2;
[x,u] = BVP1D(linspace(0,L,N),EA,0,0,0,0,0,[F,0]);

figure(1); plot(x,u)
           xlabel('position x'); ylabel('displacement u')

[u2,strain] = pwquadinterp(x,u,x);               %% evaluation at nodes
x_fine = linspace(0,L,501);
[u_fine,strain_fine] = pwquadinterp(x,u,x_fine); %% interpolation to a fine grid

figure(2); plot(x,strain,x_fine,strain_fine)
           xlabel('position x'); ylabel('strain du/dx')
           legend('at nodes', 'fine grid')

%%figure(1); print -dpdfcrop BeamStretch.pdf
%%figure(2); print -dpdfcrop BeamStretchStrain.pdf
