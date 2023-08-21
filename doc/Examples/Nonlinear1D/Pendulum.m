N = 1001; T = 2.0;
interval = linspace(0,T,N)';
BCleft = 0; BCright = [0,0];
f = {@(t,u)sin(u), @(t,u)cos(u)};
[t,u] = BVP1DNL(interval,1,0,0,1,f,BCleft,BCright,@(t)t,
                'Display','off','Tol',1e-8);

figure(1); plot(t,u); xlabel('time t'); ylabel('angle u')
v = FEM1DEvaluateDu(t,u);
disp(sprintf('For T = %g: initial angular velocity v(0) = %g, maximal angle u(T) = %g',T,v(1),u(end)))
KineticEnergy = v(1)^2/2; Potential = 1-cos(u(end));
disp(sprintf('Kinetic energy at t=0: %g, potential energy at t=T: %g, difference: %g',
              KineticEnergy, Potential, KineticEnergy-Potential))

%% determine the travel time T by an integral
u_max = u(end);
dtdu = @(u)1/sqrt(2*(cos(u)-cos(u_max)));
TravelTime = quad(dtdu,0,u_max)
