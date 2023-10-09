clear *
%% script file to solve Black-Scholes for a European put option
K = 120;          % strike
r = 0.076;        % annual gain of stock
sigma = 0.13;     % volatility
r0 = r+sigma^2/2; % safe interest rate
T = 1.0;          % maximal time to maturity
%%%%%%%%%%%%%%%
a = -0.5; b = +0.5;
interval = log(K)+linspace(a,b,101)';
BCright = 0; BCleft = K-exp(min(interval));

u0 = @(z)max(0,K-exp(z));
[z,V,tau] = IBVP1D(interval,1,sigma^2/2,-r,+r0,0,0,BCleft,BCright,u0,0,T,[10,10]);

figure(1); mesh(tau,exp(z),V); xlabel('\tau'); ylabel('S'); zlabel('V(S,\tau)')
           xlim([0,T]); ylim([100,135]); zlim([0,30]); caxis([0,25]); view([140,30])
S = linspace(100,140,101)';
V0 = interp1(exp(z),V(:,1),S); Vend = interp1(exp(z),V(:,end),S);
t_ind = find(abs(tau-T/2)<100*eps); Vmid = interp1(exp(z),V(:,t_ind),S);
figure(2); plot(S,V0,S,Vmid,S,Vend); xlabel('S'); ylabel('V(S,\tau)')
           legend('at \tau=0','at \tau=T/2','at \tau=T','location','northeast')

if 1 %% to determine the resuts with the financial package
  pkg load financial
  S = 130; tau = T/2;
  [Call,Put] = blsprice(S,K,r+sigma^2/2,tau,sigma ,0);
  disp(sprintf("For S = %g obtain Call = %g or put = %g at time tau = %g",S,Call,Put,tau))
endif

