N = 10; x = linspace(0,3,N);
[xn,u] = BVP1D(x,1,0,0,1,@(x)(1-x).^2,2,[-2,0]);
figure(1); plot(xn,u); xlabel('x'); ylabel('u')
