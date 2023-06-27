## Copyright (C) 2023 Andreas Stahel
##
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.


## Author: Andreas Stahel <andreas.stahel@gmx.com>
## Created: 2023-06-24


## -*- texinfo -*-
## @deftypefn{function file}{}[@var{x},@var{u},@var{t}] = I2BVP1D(@var{interval},@var{w1},@var{a},@var{b},@var{c},@var{d},@var{f},@var{BCleft},@var{BCright},@var{u0},@var{u1},@var{t0},@var{tend},@var{steps})
##
##   solve a second order 1D initial boundary value problem (IBVP)
##
##    w2(x)*d^2/dt^2 u(x,t) + w1(x)*d/dt u(x,t) - (a(x)*u'(x,t))' + b(x)*u'(x,t) + c(x)*u(x,t) + d(x)*f(x,t) = 0
##
##      with initial condition u(x,t0) = u0(x) and d/dt u(x,t0) = u1
##      and boundary conditions at the two endpoints
##@itemize
##@item Dirichlet: u(x,t) = g_D
##@item Neumann: a(x)*u'(x,t) = g_N1 + g_N2*u(x)
##@end itemize
##
##parameters:
##@itemize
##@item @var{interval} the discretized interval for the BVP
##@item @var{w2} constant, vector or function handle to evaluate w2(x)
##@item @var{w1} constant, vector or function handle to evaluate w1(x)
##@item @var{a} constant, vector or function handle to evaluate a(x)
##@item @var{b} constant, vector or function handle to evaluate b(x)
##@item @var{c} constant, vector or function handle to evaluate c(x)
##@item @var{d} constant, vector or function handle to evaluate d(x)
##@item @var{f} constant, vector or function handle to evaluate the f(x,t)
##@item @var{BCleft} and @var{BCright} the two boundary conditions
##@itemize
##@item for a Dirichlet condition specify a single value @var{g_D}
##@item for a Neumann condition specify the values @var{[g_N1,g_N2]}
##@end itemize
##@item @var{u0} constant, vector with the initial values at the nodes or a function handle to evaluate u(t0)
##@item @var{u1} constant, vector with the initial velocitiesat the nodes or a function handle to evaluate u(t0)
##@item @var{t0}, @var{tend} are the initial and final times
##@item @var{steps} is a vector with one or two positive integers.
##@itemize
##@item If @var{steps} = n, then n steps are taken and the n+1 results returned.
##@item If @var{steps} = [n,nint], then n*nint steps are taken and (n+1) results returned.
##@end itemize
##@end itemize
##
##return values
##@itemize
##@item @var{x} the nodes in the given interval
##@item @var{u} is a matrix with n+1 columns with the values of the solution at the nodes at different times @var{t}
##@item @var{t} is the vector with the values of the times at which the solutions are returned.
##@end itemize
##
## @end deftypefn

function  [x,u,t] = I2BVP1D(interval,w2,w1,a,b,c,d,f,BCleft,BCright,u0,u1,t0,tend,steps,varargin)

solver = 'IMPLICIT';  %% default value is IMPLICIT

if (~isempty(varargin))
  for cc = 1:2:length(varargin)
    switch tolower(varargin{cc})
      case {'solver'}
	solver = toupper(varargin{cc+1});
      otherwise
	error('Invalid optional argument, %s',varargin{cc}.name);
    endswitch % switch
  endfor % for
endif % if

[A,M,x]   = GenerateFEM1D(interval,a,b,c,d);
[x,W1,W2] = GenerateWeight1D(interval,w1,w2);

if (length(u0)==1)&&isnumeric(u0)  %% u0 given as vector
  u0 = u0*ones(size(x));
elseif ~isnumeric(u0)              %% u0 given as handle
  u0 = u0(x);
endif
if (length(u1)==1)&&isnumeric(u1)  %% u1 given as vector
  u1 = u1*ones(size(x));
elseif ~isnumeric(u1)              %% u1 given as handle
  u1 = u1(x);
endif

%% determine the boundary conditions
if length(BCleft)*length(BCright)==1  %% DD: Dirichlet at both ends
  BC = 'DD';
elseif (length(BCleft)>1)&&(length(BCright)==1) %% ND: Neumann on the left, Dirichlet on the right
  BC = 'ND';
elseif (length(BCleft)==1)&&(length(BCright)>1) %% DN: Dirichlet on the left, Neumann on the right
  BC = 'DN';
else  %% NN: Neumann on both endpoints
  BC = 'NN';
endif

switch BC
  case 'DD';
    a_left  = A(2:end-1,1);   %% first column
    a_right = A(2:end-1,end); %% last column
    A = A(2:end-1,2:end-1);
    W2 = W2(2:end-1,2:end-1); W1 = W1(2:end-1,2:end-1); M = M(2:end-1,:);
    uB = A\(-BCleft*a_left - BCright*a_right);
    u0 = u0(2:end-1); u1 = u1(2:end-1);
  case 'ND';
    a_right = A(1:end-1,end); %% last column
    A = A(1:end-1,1:end-1);
    W2 = W2(1:end-1,1:end-1); W1 = W1(1:end-1,1:end-1); M = M(1:end-1,:);
    A(1,1) += BCleft(2);
    RHS = - BCright*a_right; RHS(1) -= BCleft(1);
    uB = A\RHS;
    u0 = u0(1:end-1); u1 = u1(1:end-1);
  case 'DN';
    a_left = A(2:end,1); %% first column
    A = A(2:end,2:end);
    W2 = W2(2:end,2:end); W1 = W1(2:end,2:end); M = M(2:end,:);
    A(end,end) -= BCright(2);
    RHS = - BCleft*a_left; RHS(end) += BCright(1);
    uB = A\RHS;
    u0 = u0(2:end); u1 = u1(2:end);
  case 'NN';
    if (BCleft(1)==0)&&(BCright(1)==0)  %% only solve for the homogeneous BC if necessary
      uB = zeros(size(A,1),1);
    else   ## solve the system
      A(1,1) += BCleft(2); A(end,end) -= BCright(2);
      RHS = zeros(size(x)); RHS(1) -= BCleft(1); RHS(end) += BCright(1);
      uB = A\RHS;
    endif
endswitch
%%endif
uB = full(uB);  %% solution with f=0 and the nonzero BC

f_var = 1; %% f depends on time
if (length(f)>1)&&isnumeric(f)  %% f given as vector
  f_values = f;
  f_var = 0;  %% f does not depend on time
elseif isnumeric(f)           %% f given as scalar
  f_values = f*(ones(size(x)));
  f_var = 0;  %% f does not depend on time
endif

n_1 = steps(1);
if length(steps)>1
  n_2 = steps(2);
else
  n_2 = 1;
endif
dt = (tend-t0)/(n_1*n_2);
t = linspace(t0,tend,n_1+1);  tt = t0+dt;
u = zeros(length(uB),n_1+1);  u(:,1) = u0;

switch solver
  case 'IMPLICIT'
    Mleft    = W2+dt*W1+dt^2/4*A;  %% matrix for u(t+dt)
    Mmiddle  = 2*W2 - dt^2/2*A;   %% matrix for u(t)
    Mright   = -Mleft + 2*dt*W1;  %% matrix for u(t-dt)
    [L,U,P,Q] = lu(Mleft);  %% P*A*Q = L*U

    u_curr = u0-uB;  %% starting value
    if f_var
      f_values = f(x,t0);
    endif %% f_dep_t
    %% value at t0+Delta t
    u_new = W2\((W2-dt*W1)*dt*u1 + W2*u_curr + dt^2/2*(M*f_values-A*u_curr));

    for ii_t = 1:n_1
      for ii_2 = 1:n_2
	if f_var
	  f_values = f(x,tt);
	endif %% f_dep_t
	u_temp = Q*(U\(L\(P*(Mmiddle*u_new + Mright*u_curr + dt^2*(M*f_values)))));
	u_curr = u_new;
	u_new = u_temp;
	tt += dt;
      endfor % ii_2
      u(:,ii_t+1) = u_new + uB;
    endfor
  otherwise
    error('Invalid optional argument, %s, valid are implicit',solver);
endswitch %% solver

switch BC  %% add the Dirichlet values at the endpoints
  case 'DD'
    u = [BCleft*ones(1,n_1+1);u;BCright*ones(1,n_1+1)];
  case 'DN'
    u = [BCleft*ones(1,n_1+1);u];
  case 'ND'
    u = [u;BCright*ones(1,n_1+1)];
endswitch %% BC
