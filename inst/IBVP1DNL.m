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
## Created: 2023-11-20


## -*- texinfo -*-
## @deftypefn{function file}{}[@var{x},@var{u},@var{t}] = IBVP1DNL(@var{interval},@var{w},@var{a},@var{b},@var{c},@var{d},@var{f},@var{BCleft},@var{BCright},@var{u0},@var{t0},@var{tend},@var{steps},@var{options})
##
##   solve a 1D nonlinear initial boundary value problem (IBVP)
##
##    w(x)*d/dt u(x,t) - (a(x)*u'(x,t))' + b(x)*u'(x,t) + c(x)*u(x,t) = d(x)*f(x,t,u(x,t))
##
##      with initial condition u(x,t0) = u0(x)
##      and boundary conditions at the two endpoints
##@itemize
##@item Dirichlet: u(x,t) = g_D
##@item Neumann: a(x)*u'(x,t) = g_N1 + g_N2*u(x)
##@end itemize
##
##parameters:
##@itemize
##@item @var{interval} the discretized interval for the BVP
##@item @var{w} constant, vector or function handle to evaluate w(x)
##@item @var{a} constant, vector or function handle to evaluate a(x)
##@item @var{b} constant, vector or function handle to evaluate b(x)
##@item @var{c} constant, vector or function handle to evaluate c(x)
##@item @var{d} constant, vector or function handle to evaluate d(x)
##@item @var{f} a structure @var{f = @{@@(x,t,u),  @@(x,t,u)@}}. The two function handles evaluate f(x,t,u) and the partial derivative f_u(x,t,u).
##@item @var{BCleft} and @var{BCright} the two boundary conditions
##@itemize
##@item for a Dirichlet condition specify a single value @var{g_D}
##@item for a Neumann condition specify the values @var{[g_N1,g_N2]}
##@end itemize
##@item @var{u0} constant or vector with the initial values at the nodes or a function handle to evaluate u(x,t0)
##@item @var{t0}, @var{tend} are the initial and final times
##@item @var{steps} is a vector with one or two positive integers.
##@itemize
##@item If @var{steps} = n, then n steps are taken and the n+1 results returned.
##@item If @var{steps} = [n,nint], then n*nint steps are taken and (n+1) results returned.
##@end itemize
##@item @var{options} additional options, given as pairs name/value. 
##@itemize
##@item @var{"tol"} the tolerance for the iteration at each time step to stop, given as pair @var{[tolrel,tolabs]} for the relative and absolute tolerance. The iteration stops if the absolute or relative error is smaller than the specified tolerance. RMS (root means square) values are used. If only @var{tolrel} is specified @var{TolAbs=TolRel} is used. The default values are @var{tolrel = tolabs = 1e-5}.
##@item @var{"MaxIter"} the maximal number of iterations to be used. The default value is 10.
##@end itemize
##@end itemize
##
##return values
##@itemize
##@item @var{x} the nodes in the given interval
##@item @var{u} is a matrix with n+1 columns with the values of the solution at the nodes at times @var{t}
##@item @var{t} is the vector with the values of the times at which the solutions are returned.
##@end itemize
##
##@end deftypefn

function  [x,u,t] = IBVP1DNL(interval,w,a,b,c,d,f,BCleft,BCright,u0,t0,tend,steps,varargin)

tol = [1e-5, 1e-5];  %% default tolerance
MaxIter = 10;        %% default maximal number of iterations
Display = 'off';
if (~isempty(varargin))
  for cc = 1:2:length(varargin)
    switch tolower(varargin{cc})
      case {'tol'}
	tol = varargin{cc+1};
	if length(tol)==1
	  tol = tol*[1 1];
	endif
      case {'maxiter'}
	MaxIter = varargin{cc+1};
      case {'display'}
	Display = varargin{cc+1};
      otherwise
	error('Invalid optional argument, %s',varargin{cc}.name);
    endswitch % switch
  endfor % for
endif % if ĩsempty


[A,M,x] = GenerateFEM1D(interval,a,b,c,d);
[x,W]   = GenerateWeight1D(interval,w);

if (length(u0)==1)&&isnumeric(u0)  %% u0 given as scalar
  u0 = u0*ones(size(x));
elseif ~isnumeric(u0)              %% u0 given as handle
  u0 = u0(x);
endif

%% determine the boundary conditions
if length(BCleft)*length(BCright)==1            %% DD: Dirichlet at both ends
  BC = 'DD';
elseif (length(BCleft)>1)&&(length(BCright)==1) %% ND: Neumann on the left, Dirichlet on the right
  BC = 'ND';
elseif (length(BCleft)==1)&&(length(BCright)>1) %% DN: Dirichlet on the left, Neumann on the right
  BC = 'DN';
else                                         %% NN: Neumann on both endpoints
  BC = 'NN';
endif

f_values = f{1}(t0,x,u0);

switch BC
  case 'DD';
    a_left  = A(2:end-1,1);   %% first column
    a_right = A(2:end-1,end); %% last column
    A = A(2:end-1,2:end-1); W = W(2:end-1,2:end-1);
    M = M(2:end-1,:); Mr = M(:,2:end-1);
    uB = full(A\(-BCleft*a_left - BCright*a_right));
    %%    uB = [BCleft;uB;BCright];
    u0 = u0(2:end-1);
  case 'ND';
    a_right = A(1:end-1,end); %% last column
    A = A(1:end-1,1:end-1); W = W(1:end-1,1:end-1);
    M = M(1:end-1,:); Mr = M(:,1:end-1);
    A(1,1) += BCleft(2);
    RHS = - BCright*a_right; RHS(1) -= BCleft(1);
    uB = full(A\RHS);
    %%    uB = [uB;BCright];
    u0 = u0(1:end-1);
  case 'DN';
    a_left = A(2:end,1); %% first column
    A = A(2:end,2:end); W = W(2:end,2:end);
    M = M(2:end,:); Mr = M(:,2:end);
    A(end,end) -= BCright(2);
    RHS = - BCleft*a_left; RHS(end) += BCright(1);
    uB = full(A\RHS);
    %%    uB = [BCleft;uB];
    u0 = u0(2:end);
  case 'NN';
     Mr = M;
     if (BCleft(1)==0)&&(BCright(1)==0)  %% only solve for the homogeneous BC if necessary
       uB = zeros(size(A,1),1);
     else   ## solve the system
       A(1,1) += BCleft(2); A(end,end) -= BCright(2);
       RHS = zeros(size(x)); RHS(1) -= BCleft(1); RHS(end) += BCright(1);
       uB = full(A\RHS);
     endif
endswitch

n_1 = steps(1);
if length(steps)>1
  n_2 = steps(2);
else
  n_2 = 1;
endif
dt = (tend-t0)/(n_1*n_2);
t = zeros(n_1+1,1); t(1) = t0;  tt = t0;
u = zeros(length(uB),n_1+1);  u(:,1) = u0;

%%% CN stepper
Mleft = W + dt/2*A;
ut = u0-uB;  %% starting value
for jj1 = 1:n_1;
  for jj2 = 1:n_2
    u_new = ut;  %% starting value for the Newton iteration
    G = W*ut-0.5*dt*(A*ut - M*f_values);
    iter = 0;
    do
      iter++;
      %% evaluate f(x,t,u) and f_u(x,t,u)
      switch BC  %% add the Dirichlet values at the endpoints
	case 'DD'
	  utt = [BCleft;u_new+uB;BCright];
	  fu_values = f{2}(x,tt+dt,utt);
	  Mat_fu= M*diag(fu_values);
	  Mat_fu = Mat_fu(:,2:end-1);	  
	case 'DN'
	  utt = [BCleft;u_new+uB];
	  fu_values = f{2}(x,tt+dt,utt);
	  Mat_fu= M*diag(fu_values);
	  Mat_fu = Mat_fu(:,2:end);	  
	case 'ND'
	  utt = [u_new+uB;BCright];
	  fu_values = f{2}(x,tt+dt,utt);
	  Mat_fu= M*diag(fu_values);
	  Mat_fu = Mat_fu(:,1:end-1);	  
	case 'NN'
	  utt = u_new+uB; 
	  fu_values = f{2}(x,tt+dt,utt);
	  Mat_fu= M*diag(fu_values);
      endswitch
      f_values  = f{1}(x,tt+dt,utt);
      phi = (Mleft-dt/2*Mat_fu)\(G-(W+dt/2*A)*u_new + dt/2*M*f_values);
      u_new += phi;
      AbsError = norm(phi);
    until or(iter>=MaxIter,AbsError/sqrt(length(x))<tol(2),AbsError<tol(1)*norm(u_new))
    ut = u_new;
    tt += dt;
  endfor  %% jj2
  u(:,jj1+1) = ut+uB;  t(jj1+1) = tt;
endfor %% jj1
%%% end of CN stepper

switch BC  %% add the Dirichlet values at the endpoints
  case 'DD'
    u = [BCleft*ones(1,n_1+1);u;BCright*ones(1,n_1+1)];
  case 'DN'
    u = [BCleft*ones(1,n_1+1);u];
  case 'ND'
    u = [u;BCright*ones(1,n_1+1)];
endswitch %% BC
