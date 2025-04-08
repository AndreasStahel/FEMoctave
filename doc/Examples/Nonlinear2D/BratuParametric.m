## -*- texinfo -*-
## @deftypefn  {} {} BratuParametric.m
##
## This is a example file inside the `doc/Examples/Nonlinear2D/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

C_crit = 3.513830719;
N_C = 101; N = 201;
Interval= linspace(0,1,N)';
%% BVP1DNL, lower branch, up
C_list = linspace(1,C_crit,N_C);
u0 = @(x)0.8*sin(pi*x);
C = 1;
f = {@(x,u)C*exp(u),@(x,u)C*exp(u)};
[x,u] = BVP1DNL(Interval,1,0,0,1,f,0,0,u0,'display','iter');
u_max = max(u)
for C = C_list(2:end)
  disp(sprintf('parameter C = %f',C))
  f = {@(x,u)C*exp(u),@(x,u)C*exp(u)};
  [x,u] = BVP1DNL(Interval,1,0,0,1,f,0,0,u,'display','iter','maxiter',20);
  u_max = [u_max,max(u)];
endfor
u_max_low = u_max;

%% BVP1DNL, upper branch, up
u0 = @(x)4*sin(pi*x);
C = 1;
f = {@(x,u)C*exp(u),@(x,u)C*exp(u)};
[x,u] = BVP1DNL(Interval,1,0,0,1,f,0,0,u0,'display','iter');
u_max = max(u)
for C = C_list(2:end)
  disp(sprintf('parameter C = %f',C))
  f = {@(x,u)C*exp(u),@(x,u)C*exp(u)};
  [x,u] = BVP1DNL(Interval,1,0,0,1,f,0,0,u,'display','iter','maxiter',20);
  u_max = [u_max,max(u)];
endfor
u_max_upp = u_max;

%% BVP1DNL, lower branch, down
C_list_down = linspace(1,0.1,N_C);
u0 = @(x)0.8*sin(pi*x);
C = 1;
f = {@(x,u)C*exp(u),@(x,u)C*exp(u)};
[x,u] = BVP1DNL(Interval,1,0,0,1,f,0,0,u0,'display','iter');
u_max = max(u)
for C = C_list_down(2:end)
  disp(sprintf('parameter C = %f',C))
  f = {@(x,u)C*exp(u),@(x,u)C*exp(u)};
  [x,u] = BVP1DNL(Interval,1,0,0,1,f,0,0,u,'display','iter','maxiter',20);
  u_max = [u_max,max(u)];
endfor
u_max_low_down = u_max;
u0 = @(x)4*sin(pi*x);
C = 1;
f = {@(x,u)C*exp(u),@(x,u)C*exp(u)};
[x,u] = BVP1DNL(Interval,1,0,0,1,f,0,0,u0,'display','iter');
u_max = max(u)
for C = C_list_down(2:end)
  disp(sprintf('parameter C = %f',C))
  f = {@(x,u)C*exp(u),@(x,u)C*exp(u)};
  [x,u,inform] = BVP1DNL(Interval,1,0,0,1,f,0,0,u,'display','none','maxiter',20);
  if inform.info ==1
    u_max = [u_max,max(u)];
  else
    u_max = [u_max,NaN];
  endif
endfor
u_max_upp_down = u_max;

%% plot all
figure(2); plot(C_list_down,[u_max_low_down',u_max_upp_down'],'k',
                C_list,     [u_max_low',u_max_upp'],'k');
           xlabel('parameter C'); ylabel('u_\infty')

