## Copyright (C) 2025 Andreas Stahel
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
## @deftypefn{function file}{}[@var{u},@var{inform}] = BVP2DNL(@var{mesh},@var{a},@var{b0},@var{bx},@var{by},@var{f},@var{gD},@var{GN1},@var{gN2},@var{u0},@var{options})
##
##   solve a semilinear 2D boundary value problem
##
##@verbatim
##  -div(a*grad u - u*(bx,by))+ b0*u = f(xy,u)   in domain
##                                 u = gD        on Dirichlet boundary
##          n*(a*grad u - u*(bx,by)) = gN1+gN2*u on Neumann boundary
##@end verbatim
##
##
##parameters:
##@itemize
##@item @var{mesh} is the mesh describing the domain and the boundary types
##@item @var{a},@var{b0},@var{bx},@var{by},@var{gD},@var{gN1},@var{gN2}
##are the coefficients and functions describing the PDE.
##@*Any constant function can be given by its scalar value.
##@*The functions @var{a},@var{b0},@var{bx} and @var{by} and may also be given as vectors
##with the values of the function at the Gauss points.
##@item @var{f} function handles to evaluate f(xy) or f(xy,u) and the partial derivative with respect to u.
##@itemize
##@item @var{f = @@(xy)} a function handle to evaluate f(xy) at the Gauss points.
##@item @var{f = @{@@(xy,u),  @@(xy,u)@}} assumes that the function f depends on x,y and u. The two function handles evaluate f(xy,u) and the partial derivative f_u(xy,u).
##@end itemize
##@item @var{options} additional options, given as pairs name/value.
##@itemize
##@item @var{"tol"} the tolerance for the iteration to stop, given as pair @var{[tolrel,tolabs]} for the relative and absolute tolerance. The iteration stops if the absolute or relative error is smaller than the specified tolerance. RMS (root means square) values are used. If only @var{tolrel} is specified @var{TolAbs = TolRel} is used. The default values are @var{tolrel = tolabs = 1e-5}.
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

function  [u,inform] = BVP2DNL(Mesh,a,b0,bx,by,f,gD,gN1,gN2,u0,varargin)

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

xyGP = Mesh.GP;
nGP = size(xyGP,1);  %% number of Gauss Points

function fun_values = convert2values(fun)  %% evaluate at Gauss points
  if (length(fun)>1)&&isnumeric(fun)       %% a given as vector
    fun_values = fun;
  elseif isnumeric(fun)                    %% a given as scalar
    fun_values = fun*(ones(nGP,1));
  else                                     %% a given as function handle
    fun_values = fun(xyGP);
  endif
endfunction
a_values  = convert2values(a);
b0_values = convert2values(b0);
bx_values = convert2values(bx);
by_values = convert2values(by);

if is_function_handle(u0)
  u_values = u0(Mesh.nodes); u_valuesGauss = u0(xyGP);
elseif (length(u0)==1)&&isnumeric(u0)
  u_values = u0*ones(size(Mesh.nodes,1),1); u_valuesGauss = u0*ones(nGP,1);
else
  u_values = u0; u_valuesGauss = FEMEvaluateGP(Mesh,u_values);
endif

if isnumeric(f)  %% constant or vector of values at the nodes
  if length(f)==1
    f_values = f*ones(nGP,1);
  else
    f_values = f;
  endif%% length(f)==1
else  %% f is a function handle or cell array of function handles
  if length(f)==1  %% f(x)
    if iscell(f)
      f_values = f{1}(xyGP);
    else
      f_values = f(xyGP);
    endif %% iscell(f)
  else
    f_NL = 1;  %% f(x,u)
    if nargin(f{1})==2
      f_values = f{1}(xyGP,u_valuesGauss);
%    else
%      du = FEM1DEvaluateDu(x,u_values);
%      f_values = f{1}(x,u_values,du);
    endif %% nargin(f{1})
  endif %% length(f)==1
endif


%% solve a system with the actual BC
switch Mesh.type
  case 'linear' %% first order elements
    [A,b] = FEMEquation    (Mesh,a_values,b0_values,bx_values,by_values,f_values,gD,gN1,gN2);
    [A0,gD0] = FEMEquation(Mesh,a_values,b0_values,bx_values,by_values,0,gD,0,0);
  case 'quadratic'  %% second order elements
    [A,b] = FEMEquationQuadM(Mesh,a_values,b0_values,bx_values,by_values,f_values,gD,gN1,gN2);
    [A0,gD0] = FEMEquationQuadM(Mesh,a_values,b0_values,bx_values,by_values,0,gD,gN1,gN2);
  case 'cubic'  %% third order elements
    [A,b] = FEMEquationCubic(Mesh,a_values,b0_values,bx_values,by_values,f_values,gD,gN1,gN2);
    [A0,gD0] = FEMEquationCubic(Mesh,a_values,b0_values,bx_values,by_values,0,gD,gN1,gN2);
endswitch
u_values = FEMSolve(Mesh,A,b,gD);  %% solve the resulting system

W = GenerateWeight2D(Mesh);

relError = 2*tol(1); AbsError = 2*tol(2);
inform.info = 1; inform.iter = 0;

function phi = SolvePhi(b0,f)
  %% solve a system with the actual BC
 switch Mesh.type
   case 'linear' %% first order elements
     A = FEMEquation    (Mesh,a_values,b0,bx_values,by_values,0,0,0,gN2);
   case 'quadratic'  %% second order elements
     A = FEMEquationQuad(Mesh,a_values,b0,bx_values,by_values,0,0,0,gN2);
   case 'cubic'  %% third order elements
     A = FEMEquationCubic(Mesh,a_values,b0,bx_values,by_values,0,0,0,gN2);
 endswitch
 phi = FEMSolve(Mesh,A,f+gD0,0);  %% solve the resulting system
endfunction

ind_free = find(Mesh.node2DOF>0);  %% indices of DOF
sqrt_n = sqrt(size(Mesh.nodes,1));

if (iscell(f)&(length(f)>1))  %% only iterate for nonlinear functions
do   %% until error small enough, or inform.iter > MaxIter
  u_old = u_values;  %% solution at nodes
  inform.iter++;
  RHS = A0*u_old(ind_free);
  u_valuesGauss = FEMEvaluateGP(Mesh,u_values);
  f_values         = f{1}(xyGP,u_valuesGauss);  %% f(x,u) at GP
  dfdu_valuesGauss = f{2}(xyGP,u_valuesGauss);  %% df/du(x,u) at GGP
  phi = SolvePhi(b0_values-dfdu_valuesGauss,RHS-W*f_values);
  u_values += phi;
  AbsError = norm(u_values-u_old);
  if strcmp(Display,'iter')
%    disp(sprintf('iteration %i, RMS(correction) = %e, RMS(phi) = %e',inform.iter,AbsError/sqrt_n,norm(phi)/sqrt_n))
    disp(sprintf('iteration %i, RMS(correction) = %e',inform.iter,AbsError/sqrt_n))
  endif  %% strcmp

until or(inform.iter>=MaxIter,AbsError/sqrt_n<tol(2),AbsError<tol(1)*norm(u_values))
endif %%nargin(f{1})>1

u = u_values;  %% assign the return value u

inform.AbsError = AbsError/sqrt_n;
if (inform.iter >= MaxIter)
  disp('BVP2DNL did not converge')
  inform.info = -1;
endif
endfunction
