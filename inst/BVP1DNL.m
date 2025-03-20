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
## Created: 2023-08-21

## -*- texinfo -*-
## @deftypefn{function file}{}[@var{x},@var{u},@var{inform}] = BVP1DNL(@var{interval},@var{a},@var{b},@var{c},@var{d},@var{f},@var{BCleft},@var{BCright},@var{u0},@var{options})
##
##   solve a nonlinear 1D boundary value problem (BVP)
##
##      -(a(x,u,u')*u'(x))' + b(x)*u'(x) + c(x)*u(x) = d(x)*f(x,u,u')
##
##      with boundary conditions at the two endpoints
##@itemize
##@item Dirichlet: u(x) = g_D
##@item Neumann: a(x,u,u')*u'(x) = g_N1 + g_N2*u(x)
##@end itemize
##
##parameters:
##@itemize
##@item @var{interval} the discretized interval for the BVP
##@item @var{a} constant, vector or function handle to evaluate a(x), a(x,u) or a(x,u,u') at the Gauss points.
##@itemize
##@item @var{a} a constant or vector of values of a(x).
##@item @var{a = @@(x)} a function handle to evaluate f(x).
##@item @var{a = @@(x,u)} assumes that the function a(x,u) depends on x and u. 
##@item @var{a = @@(x,u,u')} assumes that a(x,u,u') depends on x,  u and u'.
##@end itemize
##@item @var{b} constant, vector or function handle to evaluate b(x) at Gauss points
##@item @var{c} constant, vector or function handle to evaluate c(x) at Gauss points
##@item @var{d} constant, vector or function handle to evaluate d(x) at Gauss points
##@item @var{f} constant, vector or function handle to evaluate f(x), f(x,u) or f(x,u,u') and the partial derivatives at nodes
##@itemize
##@item @var{f} a constant or vector of values of f(x) at the nodes.
##@item @var{f = @@(x)} a function handle to evaluate f(x) at the nodes.
##@item @var{f = @{@@(x,u),  @@(x,u)@}} assumes that the function f depends on x and u. The two function handles evaluate f(x,u) and the partial derivative f_u(x,u).
##@item @var{f = @{@@(x,u,u'),  @@(x,u,u') , @@(x,u,u')@}} assumes that f depends on x,  u and u'. The three function handles evaluate f(x,u,u') and the partial derivatives f_u(x,u,u') and f_u'(x,u,u').
##@end itemize
##@item @var{BCleft} and @var{BCright} the two boundary conditions
##@itemize
##@item for a Dirichlet condition specify a single value @var{g_D}
##@item for a Neumann condition specify the two values @var{[g_N1,g_N2]}
##@end itemize
##@item @var{u0} constant, vector or function handle to evaluate u0(x) at the nodes. This is the starting value for the iteration.
##@item @var{options} additional options, given as pairs name/value.
##@itemize
##@item @var{"tol"} the tolerance for the iteration to stop, given as pair @var{[tolrel,tolabs]} for the relative and absolute tolerance. The iteration stops if the absolute or relative error is smaller than the specified tolerance. RMS (root means square) values are used. If only @var{tolrel} is specified @var{TolAbs=TolRel} is used. The default values are @var{tolrel = tolabs = 1e-5}.
##@item @var{"MaxIter"} the maximal number of iterations to be used. The default value is 10.
##@item @var{"Display"} should information be displayed for the iterations
##@itemize
##@item @var{"off"} no display, default
##@item @var{"iter"} display the number of the iteration and the RMS size of the update
##@end itemize
##@end itemize
##@end itemize
##
##return values
##@itemize
##@item @var{x} the nodes in the given interval
##@item @var{u} the values of the solution at the nodes
##@item @var{inform} a structure with information on the performance of the algorithm
##@itemize
##@item @var{inform.info} = 1 if the algorithm converged with the desired tolerance, -1 if not.
##@item @var{inform.iter} the number of iterations used.
##@item @var{inform.AbsError} the RMS value of the last correction applied.
##@end itemize
##@end itemize
##
##@end deftypefn

function  [x,u,inform] = BVP1DNL(interval,a,b,c,d,f,BCleft,BCright,u0,varargin)

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

x = [interval(:) interval(:)+[diff(interval(:))/2;0]]'(1:end-1); x = x(:);
[xGauss, Nodes2GaussU, Nodes2GaussDU] = FEM1DGaussPoints(x);
n = length(x);  %% number of nodes

function fun_values = convert2values(fun)  %% evaluate at Gauss points
  if (length(fun)>1)&&isnumeric(fun)       %% a given as vector
    fun_values = fun;
  elseif isnumeric(fun)                    %% a given as scalar
    fun_values = fun*(ones(size(xGauss)));
  else                                     %% a given as function handle
    fun_values = fun(xGauss);
  endif
endfunction
b_values = convert2values(b);
c_values = convert2values(c);
d_values = convert2values(d);

if is_function_handle(u0)
  u_values = u0(x); u_valuesGauss = u0(xGauss);
elseif (length(u0)==1)&&isnumeric(u0)
  u_values = u0*ones(size(x)); u_valuesGauss = u0*ones(size(xGauss));
else
  u_values = u0; u_valuesGauss = Nodes2GaussU*u_values;
endif

a_NL = 0;  %% flag if coefficient a depends on u and u'
if isnumeric(a)
  a_values = convert2values(a);
else
  switch nargin(a)
    case 1  %% a(x) depends on x only
      a_values = convert2values(a);
    case 2  %% a(x,u) depends on x and u
      a_values = a(xGauss,u_valuesGauss);
      a_NL = 1;
    case 3 %% a(x,u,u') depends on x, u and u'
      du_valuesGauss = Nodes2GaussDU*u_values;
      a_values = a(xGauss,u_valuesGauss,du_valuesGauss);
      a_NL = 1;
  endswitch
endif

f_NL = 0;  %% flag if function f depends on u and u'
if isnumeric(f)  %% constant or vector of values at the nodes
  if length(f)==1
    f_values = f*ones(size(x));
  else
    f_values = f;
  endif%% length(f)==1
else  %% f is a function handle or cell array of function handles
  if length(f)==1  %% f(x)
    if iscell(f)
      f_values = f{1}(x);
    else
      f_values = f(x);
    endif %% iscell(f)
  else
    f_NL = 1;  %% f(x,u) or f(x,u,u')
    if nargin(f{1})==2
      f_values = f{1}(x,u_values);
    else
      du = FEM1DEvaluateDu(x,u_values);
      f_values = f{1}(x,u_values,du);
    endif %% nargin(F{1})
  endif %% length(f)==1
endif


%% determine the boundary conditions
if length(BCleft)*length(BCright)==1  %% DD: Dirichlet at both ends
  BC = "DD";
elseif (length(BCleft)>1)&&(length(BCright)==1) %% ND: Neumann on the left, Dirichlet on the right
  BC = "ND";
elseif (length(BCleft)==1)&&(length(BCright)>1) %% DN: Dirichlet on the left, Neumann on the right
  BC = "DN";
else  %% NN: Neumann on both endpoints
  BC = "NN";
endif

%% apply the BC to the matrix
function [a_left,a_right] = BCcontribA(A)
  switch BC
    case "DD"
      a_left  = A(2:end-1,1);   %% first column
      a_right = A(2:end-1,end); %% last column
    case "ND";
      a_right = A(1:end-1,end); %% last column
      a_left = 0;
    case "DN";
      a_left = A(2:end,1); %% first column
      a_right = 0;
    case"NN";
      a_left = 0; a_right = 0;
  endswitch
endfunction

function An = ApplyBC2MatrixA(A)
  switch BC
    case "DD"
      An = A(2:end-1,2:end-1);
    case "ND";
      An = A(1:end-1,1:end-1);
      An(1,1) += BCleft(2);
    case "DN";
      An = A(2:end,2:end);
      An(end,end) -= BCright(2);
    case"NN";
      An(1,1) += BCleft(2); An(end,end) -= BCright(2);
  endswitch
endfunction

function M = ApplyBC2MatrixM(M)
  switch BC
    case "DD"
      M = M(2:end-1,:);
    case "ND";
      M = M(1:end-1,:);
    case "DN";
      M = M(2:end,:);
  endswitch
endfunction

%% solve a system with the actual BC
function sol = SolveSystem(A,RHS,BCleft,BCright,a_left,a_right)
  switch BC
  case "DD";
    sol = A\(RHS - BCleft*a_left - BCright*a_right);
    sol = [BCleft;sol;BCright];
  case "ND";
    RHS = RHS - BCright*a_right; RHS(1) -= BCleft(1);
    sol = [A\RHS;BCright];
  case "DN";
    RHS = RHS - BCleft*a_left; RHS(end) += BCright(1);
    sol = [BCleft;A\RHS];
  case "NN";
    RHS(1) -= BCleft(1); RHS(end) += BCright(1);
    sol = A\RHS;
  endswitch
endfunction

%% solve a system with homogeneous BC
function sol = SolveSystemPhi(A,RHS,BCleft,BCright)
  switch BC
  case "DD";
    sol = A\(RHS);
    sol = [0;sol;0];
  case "ND";
    RHS(1) -= BCleft(1);
    sol = [A\RHS;0];
  case "DN";
    RHS(end) += BCright(1);
    sol = [0;A\RHS];
  case "NN";
    RHS(1) -= BCleft(1); RHS(end) += BCright(1);
    sol = A\RHS;
  endswitch
endfunction

%%% solve the linear BVP once, assure the correct BC
[A,M,x] = GenerateFEM1D(interval,a_values,b_values,c_values,d_values);

A0 = A;   %% no BC applied
A00 = ApplyBC2MatrixA(A0);
[a_left,a_right] = BCcontribA(A);
A = ApplyBC2MatrixA(A); M = ApplyBC2MatrixM(M);
u_values  = SolveSystem(A,M*f_values,BCleft,BCright,a_left,a_right);

relError = 2*tol(1); AbsError = 2*tol(2);
inform.info = 1; inform.iter = 0;

do   %% until error small enough, or inform.iter > MaxIter
  u_old = u_values;
  inform.iter++;
  if inform.iter == 1  %% initiate A2
    %%A2 = A0;
    A2 = ApplyBC2MatrixA(A0);
  endif %% inform.iter

  u_valuesGauss  = Nodes2GaussU  * u_values;
  du_valuesGauss = Nodes2GaussDU * u_values;

  if a_NL ;  %% update the coefficient a, if necessary
    switch nargin(a)
      case 2  %% a depends on x and u
	a_values = a(xGauss,u_valuesGauss);
      case 3 %% a depends on x, u and u'
	a_values = a(xGauss,u_valuesGauss,du_valuesGauss);
    endswitch
  endif  % a_NL

  if f_NL   %%% evaluate f and apply one Newton step
    switch length(f)  %% evaluate f and its partial derivatives
      case 1   %% f(x)
        %%f_values = f(x);                              %% f(x) at nodes
        A2 = GenerateFEM1D(interval,a_values,b_values,c_values,d_values);
      case 2   %% f(x,u)
	f_values         = f{1}(x,u_values);            %% f(x,u) at nodes
	dfdu_valuesGauss = f{2}(xGauss,u_valuesGauss);  %% df/du(x,u) at Gauss
	A2 = GenerateFEM1D(interval,a_values,b,c_values-d_values.*dfdu_valuesGauss,d);
      case 3   %% f(x,u,u')
	du_values = FEM1DEvaluateDu(x,u_values);     %% u' at nodes
	f_values  = f{1}(x,u_values,du_values);      %% f(x,u,u') at nodes
	dfdu_valuesGauss   = f{2}(xGauss,u_valuesGauss,du_valuesGauss);% df/du
	dfdupr_valuesGauss = f{3}(xGauss,u_valuesGauss,du_valuesGauss);% df/du'
	A2 = GenerateFEM1D(interval,a_values,
			   b_values-d_values.*dfdupr_valuesGauss,
			   c_values-d_values.*dfdu_valuesGauss,d);
    endswitch  %% length(f)

    A2 = ApplyBC2MatrixA(A2);
    if or(inform.iter==1,a_NL==0)
      Au = A0*u_values;
    else
      Au = A*u_values;
    endif

    switch BC  %% solve the problem for phi and make one Newton step
      case "DD"
	phi = SolveSystemPhi(A2,-Au(2:end-1)+M*f_values,BCleft,BCright);
      case "ND"
	phi = SolveSystemPhi(A2,-Au(1:end-1)+M*f_values,BCleft,BCright);
      case "DN"
	phi = SolveSystemPhi(A2,-Au(2:end)+M*f_values,BCleft,BCright);
      case "NN"
	phi = SolveSystemPhi(A2,-Au+M*f_values,BCleft,BCright);
    endswitch
    u_values += phi;
    %%% evaluate f(x,u) or f(x,u,u')
    switch length(f)  %% evaluate f and its partial derivatives
      case 2   %% f(x,u)
	f_values  = f{1}(x,u_values);            %% f(x,u) at nodes
      case 3   %% f(x,u,u')
	du_values = FEM1DEvaluateDu(x,u_values); %% u' at nodes
	f_values  = f{1}(x,u_values,du_values);  %% f(x,u,u') at nodes
    endswitch  %% length(f)
  endif %% f_NL

  if a_NL %% apply successive substitution step with new matrix
    u_valuesGauss = Nodes2GaussU*u_values;
    switch nargin(a)
      case 2 %% a(x,u) depends on x and u
	a_values = a(xGauss,u_valuesGauss);
      case 3 %% a(x,u,u') depends on x, u and u'
	du_valuesGauss = Nodes2GaussDU*u_values;
	a_values = a(xGauss,u_valuesGauss,du_valuesGauss);
    endswitch
    if f_NL   %%% evaluate f(x,u) or f(x,u,u')
      switch length(f)  %% evaluate f and its partial derivatives
	case 2   %% f(x,u)
	  f_values  = f{1}(x,u_values);            %% f(x,u) at nodes
	case 3   %% f(x,u,u')
	  du_values = FEM1DEvaluateDu(x,u_values); %% u' at nodes
	  f_values  = f{1}(x,u_values,du_values);  %% f(x,u,u') at nodes
      endswitch  %% length(f)
    endif %% f_NL

    A = GenerateFEM1D(interval,a_values,b_values,c_values,d_values);

    [a_left,a_right] = BCcontribA(A);  A_BC = ApplyBC2MatrixA(A);
    u_values  = SolveSystem(A_BC,M*f_values,BCleft,BCright,a_left,a_right);

  else     %% if a_NL  %%  solve the system with the old matrix
    u_values = SolveSystem(A00,M*f_values,BCleft,BCright,a_left,a_right);
  endif %% a_NL %% apply successive substitution step

  AbsError = norm(u_values-u_old);
  if strcmp(Display,'iter')
    if f_NL
      disp(sprintf('iteration %i, RMS(correction) = %e, RMS(phi) = %e',inform.iter,AbsError/sqrt(n),norm(phi)/sqrt(n)))
    else
      disp(sprintf('iteration %i, RMS(correction) = %e',inform.iter,AbsError/sqrt(n)))
    endif %% f_NL
  endif  %% strcmp

  u = u_values;  %% assign the return value u

  if f_NL
    AbsError2 = mean([AbsError,norm(phi)]);
  else
    AbsError2 = AbsError;
  endif %% f_NL

until or((a_NL+f_NL)==0,inform.iter>=MaxIter,AbsError2/sqrt(n)<tol(2),AbsError2<tol(1)*norm(u))
inform.AbsError = AbsError/sqrt(n);
if (inform.iter >= MaxIter)
  disp('BVP1DNL did not converge')
  inform.info = -1;
endif
endfunction
