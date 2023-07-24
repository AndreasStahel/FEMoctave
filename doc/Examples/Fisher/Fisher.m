function Fisher()
c = 5/sqrt(6);
M = 20; N = 50;
f = {@(x,u)(u.*(1-u)), @(x,u)(1-2*u)};
Join = 0.25;
interval = linspace(0,M,N)';  intervaln = linspace(-M,0,N)';
%% find the parameter
function slopes = FindSlopes(valueM,graphs)
 BCleft = Join; BCright = valueM;
 [xp,u0] = BVP1D(interval,1,-c,0,1,0,BCleft,BCright);
 [xp,up,inform] = BVP1DNL(interval,1,-c,0,1,f,BCleft,BCright,u0);
 dup = FEM1DEvaluateDu(xp,up);

 BCleft = 1; BCright = Join;
 [xn,u0] = BVP1D(intervaln,1,-c,0,1,0,BCleft,BCright);
 [xn,un,inform] = BVP1DNL(intervaln,1,-c,0,1,f,BCleft,BCright,u0);
 dun = FEM1DEvaluateDu(xn,un);

 slopes = [dun(end),dup(1)];
endfunction

%%ValueMList = linspace(0,1e-3,10);  %% for M = 10
ValueMList = linspace(0,2e-7,10);  %% for M = 20
for jj = 1:length(ValueMList)
 s = FindSlopes(ValueMList(jj));
 sn(jj) = s(1); ,sp(jj) = s(2);
endfor

figure(1); plot(ValueMList,sp,ValueMList,sn);
xlabel('value of u(+M)'); ylabel('slopes  at x=0');
legend('x>0','x<0','location','south')

%% find the best value at x=+M
SameSlope = @(valueM)diff(FindSlopes(valueM));
%%ValueAtM = fzero(SameSlope,[2,3]*1e-4)   %% for M = 10
ValueAtM = fzero(SameSlope,[5,9]*1e-8)     %% for M = 20

%% solve the two BVPs
BCleft = Join; BCright = ValueAtM;     %% for M = 20
%%BCright = 2.7509e-4;  %% for M = 10
u0 = @(x)(xp-M).^6/M^6*Join;
[xp,up,inform] = BVP1DNL(interval,1,-c,0,1,f,BCleft,BCright,u0,
                         'Display','off','tol',1e-8);
inform
figure(11); plot(xp,u0(xp),xp,up); xlabel('x'); ylabel('u');
            legend('u_0','u','location','northeast')

intervaln = linspace(-M,0,N)'; BCleft = 1; BCright = Join;
u0 = @(x)1-(x+M).^4/M^4*(1-Join);

[xn,un,inform] = BVP1DNL(intervaln,1,-c,0,1,f,BCleft,BCright,u0,
                         'Display','off','tol',1e-8);
inform
figure(12); plot(xn,u0(xn),xn,un); xlabel('x'); ylabel('u');
            legend('u_0','u','location','southwest')

figure(2); plot(xn,un,xp,up); xlabel('x'); ylabel('u');

%% evaluate derivatives at x=0
[dup, ddup] = FEM1DEvaluateDu(xp,up); [dun, ddun] = FEM1DEvaluateDu(xn,un);
disp(sprintf("first derivatives: u'(x-0)= %g  and u'(x+0)= %g, difference = %g",
              dun(end),dup(1),dun(end)-dup(1)))
endfunction
