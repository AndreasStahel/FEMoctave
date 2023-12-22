%% see AgarHodiRega19 p. 277 Problem 8.13
m = 1; g = 9.81; T = 2; Interval = linspace(-1,1,51)';
f = {@(x,u,du)m*g/T*sqrt(1+du.^2), @(x,u,du)0, @(x,u,du)m*g/T*du./sqrt(1+du.^2)};

[x,u,inform] = BVP1DNL(Interval,1,0,0,-1,f,0,0,0, 'Tol',1e-12, 'Display','iter');
inform
figure(1); plot(x,u); xlabel('x'); ylabel('u_{FEM}(x)');

u_exact = @(x)T/(m*g)*(-cosh(m*g/T) + cosh(m*g/T*x));
figure(2); plot(x,u-u_exact(x)); xlabel('x'); ylabel('u_{FEM} - u_{exact}');

du = FEM1DEvaluateDu(x,u);
Length = [trapz(x,sqrt(1+du.^2)),2*T/(m*g)*sinh(m*g/T)]
