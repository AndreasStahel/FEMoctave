## Copyright (C) 2020 Andreas Stahel
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
## Created: 2020-03-30

function [u,t] = IBVP2DNL(Mesh,m,a,b0,bx,by,f,gD,gN1,gN2,u0,t0,tend,steps,varargin)
  ## -*- texinfo -*-
  ## @deftypefn{function file}{}[@var{u},@var{t}] = IBVP2D(@var{Mesh},@var{m},@var{a},@var{b0},@var{bx},@var{by},@var{f},@var{gD},@var{gN1},@var{gN2},@var{u0},@var{t0},@var{tend},@var{steps},@var{options})
  ##
  ##   Solve an initial boundary value problem
  ##
  ##@verbatim
  ## m*d/dt u - div(a*grad u-u*(bx,by)) + b0*u = f(xy,t,u) in domain
  ##                                         u = gD        on Dirichlet boundary
  ##                   n*(a*grad u -u*(bx,by)) = gN1+gN2*u on Neumann boundary
  ##                                     u(t0) = u0        initial value
  ##@end verbatim
  ##
  ##parameters:
  ##@itemize
  ##@item @var{Mesh} is the mesh describing the domain and the boundary types
  ##@item @var{m},@var{a},@var{b0},@var{bx},@var{by},@var{f},@var{gD},@var{gN1},@var{gN2}
  ##are the coefficients and functions describing the PDE.
  ##@*Any constant function can be given by its scalar value.
  ##@*The functions @var{m},@var{a},@var{b0},@var{bx} and @var{by} may also be given as vectors with the values of the function at the Gauss points.
  ##@item @var{f} may be given as a string or handle for a function depending on (x,y), time t amd u or as a vector with the values at nodes or as scalar.
  ##If @var{f} is given by a scalar or vector it is independent on time t and u.
  ##@item @var{u0} is the initial value, can be given as a constant, function name or as vector with the values at the nodes
  ##@item @var{t0}, @var{tend} are the initial and final times
  ##@item @var{steps} is a vector with one or two positive integers.
  ##@itemize
  ##@item If @var{steps} = n, then n Crank Nicolson steps are taken and the n+1 results returned.
  ##@item If @var{steps} = [n,nint], then n*nint steps are taken and (n+1) results returned.
  ##@end itemize
  ##@item @var{options} additional options, given as pairs name/value.
  ##@itemize
  ##@item The stepping algorithm can be selected as @var{"solver"}. The possible values
  ##@itemize
  ##@item @var{"CNPC"} standard Crank-Nicolson with predictor corrector (default)
  ##@item @var{"CNEXP"} standard Crank-Nicolson with explicit nonlinearity
  ##@end itemize
  ##@item Complex coefficients can be selected by @var{type}.
  ##The possible values are
  ##@itemize
  ##@item @var{"real"}: all coefficients are real (default)
  ##@item @var{"complex"}: some coefficients might be complex
  ##@end itemize
  ## If only the coefficient @var{m} is complex, there is no need to ask for complex coefficients.
  ##@end itemize
  ##@end itemize
  ##
  ##return values
  ##@itemize
  ##@item @var{u} is a matrix with n+1 columns with the values of the solution at the nodes at different times @var{t}
  ##@item @var{t} is the vector with the values of the times at which the solutions are returned.
  ##@end itemize
  ##
  ## @c Will be cut out in ??? info file and replaced with the same
  ## @c references explicitly there, since references to core Octave
  ## @c functions are not automatically transformed from here to there.
  ## @c BEGIN_CUT_TEXINFO
  ## @seealso{I2BVP2D, BVP2D, BVP2Dsym, BVP2eig, CreateMeshRect, CreateMeshTriangle}
  ## @c END_CUT_TEXINFO
  ## @end deftypefn


%  [u,t] = IBVP2D(Mesh,m,a,b0,bx,by,f,gD,gN1,gN2,u0,t0,tend,steps,options)
%  Solve an initial boundary value problem
%
%  m*d/dt u -div(a*grad u-u*(bx,by)) + b0*u = f         in domain
%                                         u = gD        on Dirichlet boundary
%                              n*(a*grad u) = gN1+gN2*u on Neumann boundary
%                                     u(t0) = u0        initial value
%
%   [u,t] = IBVP2D(Mesh,'m','a','b0','f','gD','gN1','gN2','u0',t0,tend,steps)
%   [u,t] = IBVP2D(Mesh,mVec,aVec,b0Vec,fVec,'gD','gN1','gN2',u0Vec,t0,tend,steps)
%
%      Mesh is the mesh describing the domain
%      m,a,b0,bx,by,f,gD,gN1,gN2
%          are the names of the functions or scalar values
%          The functions m, a, b, bx and by may also be given as vector
%          with the values of the function at the Gauss points
%        f may be given as a string or handle for a function depending on (x,y), time t and u
%          or a a vector with the values at nodes.
%          If f is given by a scalar or vector it is independent on time t and u
%       u0 the initial value can be given as a constant, function name
%          or as vector of the values at the nodes
%       t0 is the initial time
%     tend is the final time
%    steps is the number of Crank Nicolson steps to be taken
%          if steps=[n,m], the n*m steps of equal length will be taken,
%          and n intermediate results are returned
%  options additional options, given as pairs name/value.
%          Currently only the stepping algorithm can be selected as "solver"
%          and the possible values
%          "CNPC"     standard Crank-Nicolson with predictor-corrector (default)
%          "CNEXP"    standard Crank-Nicolson with explicity nolinearity
%
%        u is the matrix with n+1 columns, each of them containing
%          the solution at a time
%        t = linspace(t0,tend,steps(1)+1) are the times at which
%          the solutions are returned
%
% see also BVP2D, BVP2Dsym, BVP2Deig

solver = 'CNPC'; %% default value is Crank-Nicolson with predictor corrector
Type   = 'REAL'; %% default value
if (~isempty(varargin))
  for cc = 1:2:length(varargin)
    switch toupper(varargin{cc})
      case {'TYPE'}
	Type = toupper(varargin{cc+1});
      case {'SOLVER'}
	solver = toupper(varargin{cc+1});
      otherwise
	error('Invalid optional argument, %s. Possible values: SOLVER, TYPE',varargin{cc});
    endswitch % switch
  endfor % for
endif % if isempty

if ((strcmp(Type,'REAL')==0)&&(strcmp(Type,'COMPLEX')==0))
  error('wrong TYPE, possible values are REAL or COMPLEX')
endif


Mesh.node2DOF = Mesh.node2DOF(:,1);  %% not an elasticity problem

if (gD==0)&&(gN1==0)  %% only solve for the homogeneous BC if necessary
  u_B = 0;
else
  u_B = BVP2D(Mesh,a,b0,bx,by,0,gD,gN1,0,'type',TYPE);  %% solve BVP
endif
switch Mesh.type
case 'linear'    %% linear elements
  switch Type
    case 'REAL'
      if exist("FEMEquation.oct")==3
	A = FEMEquation(Mesh,a,b0,bx,by,0, 0, 0,gN2);;
      else
	A = FEMEquationM(Mesh,a,b0,bx,by,0, 0, 0,gN2);;
      endif
    case 'COMPLEX'
      if exist("FEMEquationComplex.oct")==3
	A = FEMEquationComplex(Mesh,a,b0,bx,by,0, 0, 0,gN2);;
      else
	A = FEMEquationM(Mesh,a,b0,bx,by,0, 0, 0,gN2);;
      endif
  endswitch
case 'quadratic' %% quadratic elements
  switch Type
    case 'REAL'
      if exist("FEMEquationQuad.oct")==3
	A = FEMEquationQuad(Mesh,a,b0,bx,by,0, 0, 0,gN2);;
      else
	A = FEMEquationQuadM(Mesh,a,b0,bx,by,0, 0, 0,gN2);;
      endif
    case 'COMPLEX'
      if exist("FEMEquationQuadComplex.oct")==3
	A = FEMEquationQuadComplex(Mesh,a,b0,bx,by,0, 0, 0,gN2);;
      else
	A = FEMEquationQuadM(Mesh,a,b0,bx,by,0, 0, 0,gN2);;
      endif
  endswitch
case 'cubic'     %% cubic elements
  switch Type
    case 'REAL'
      if exist("FEMEquationCubic.oct")==3
	A = FEMEquationCubic(Mesh,a,b0,bx,by,0, 0, 0,gN2);;
      else
	A = FEMEquationCubicM(Mesh,a,b0,bx,by,0, 0, 0,gN2);;
      endif
    case 'COMPLEX'
      if exist("FEMEquationCubicComplex.oct")==3
	A = FEMEquationCubicComplex(Mesh,a,b0,bx,by,0, 0, 0,gN2);;
      else
	A = FEMEquationCubicM(Mesh,a,b0,bx,by,0, 0, 0,gN2);;
      endif
  endswitch
endswitch

Wu  = FEMInterpolWeight(Mesh,m);  %% weight matrix, leading to W*(d/dt u)
Wf  = FEMInterpolWeight(Mesh,1);  %% weight matrix, leading to W*f

if length(steps)==1
  dt = (tend-t0)/steps;
  steps = [steps,1];
else
  dt = (tend-t0)/(steps(1)*steps(2));
endif

if ischar(u0)
  u0 = feval(u0,Mesh.nodes);
elseif isscalar(u0)
  u0 = u0*ones(length(Mesh.nodesT),1);
else
  u0 = u0(:);
endif

ind_free = find(Mesh.node2DOF(:,1)>0);       %% which nodes lead to DOF
ind_Dirichlet = find(Mesh.node2DOF(:,1)==0); %% Dirichlet nodes
W = Wu(:,ind_free);

t = t0;  f_dep_t = false;
if ischar(f)
  fVec = feval(f,Mesh.nodes,t+dt/2,u0);
  f_dep_t = true;  % has to be evaluated at each timestep
elseif is_function_handle(f)
  fVec = f(Mesh.nodes,t+dt/2,u0);
  f_dep_t = true;  % has to be evaluated at each timestep
elseif isscalar(f)
  fVec = f*ones(length(Mesh.nodesT),1);
else
  fVec = f;
endif

switch solver
 case 'CNEXP'
 Mleft = W+dt/2*A;  Mright = W-dt/2*A;
 [L,U,P,Q,R] = lu(Mleft);  %% P*(R\A)*Q = L*U
 t = t0;
 u = zeros(length(u0),steps(1)+1);

 u_new = u0-u_B; u_new(ind_Dirichlet) = 0;
 u(:,1) = u0;
 for ii_t = 1:steps(1)
   for ii_2 = 1:steps(2)
     if f_dep_t
       fVec = feval(f,Mesh.nodes,t+dt/2,u_new+u_B);
     endif %% f_dep_t
     u_temp = Q*(U\(L\(P*(R\(Mright*u_new(ind_free) + dt*(Wf*fVec))))));
     u_new(ind_free) = u_temp;
     t += dt;
   endfor % ii_2
   u(:,ii_t+1) = u_new + u_B;
 endfor

 case 'CNPC'
 Mleft = W+dt/2*A;  Mright = W-dt/2*A;
 [L,U,P,Q,R] = lu(Mleft);  %% P*(R\A)*Q = L*U
 t = t0;
 u = zeros(length(u0),steps(1)+1);

 u_new  = u0-u_B; u_new (ind_Dirichlet) = 0;
 u_new1 = u0-u_B; u_new1(ind_Dirichlet) = 0;
 u(:,1) = u0;
 for ii_t = 1:steps(1)
   for ii_2 = 1:steps(2)
     if f_dep_t
       fVec1 = feval(f,Mesh.nodes,t,u_new+u_B);
     endif %% f_dep_t
     u_temp = Q*(U\(L\(P*(R\(Mright*u_new(ind_free) + dt*(Wf*fVec1))))));
     u_new1(ind_free) = u_temp;
     if f_dep_t
       fVec2 = feval(f,Mesh.nodes,t,u_new1+u_B);
     endif %% f_dep_t
     u_temp = Q*(U\(L\(P*(R\(Mright*u_new(ind_free) + dt*(Wf*(fVec1+fVec2)/2))))));
     u_new(ind_free) = u_temp;
     t += dt;
   endfor % ii_2
   u(:,ii_t+1) = u_new + u_B;
 endfor
 otherwise
    error('Invalid optional argument for solver: %s, valid are CNPC, CNEXP',solver);
endswitch
t = linspace(t0,tend,steps(1)+1);
endfunction
