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
## Created: 2023-10-26


## -*- texinfo -*-
## @deftypefn{function file}{}[@var{u1},@var{u2},@var{t}] = PlaneStressDynamic(@var{mesh},@var{E},@var{nu},@var{rho},@var{f},@var{gD},@var{gN},@var{u0},@var{v0},@var{t0},@var{tend},@var{steps},@var{options})
##
##   solve a dynamic plane stress problem
##
##
##parameters:
##@itemize
##@item @var{mesh} is the mesh describing the domain and the boundary types
##@item @var{E},@var{nu} Young's modulus and Poisson's ratio for the material
##@item @var{rho},the density of the material
##@item @var{f = @{f1,f2@}} a cell array with the two components of the volume forces. The functions take two arguments, coordinates xy and time t, i.e. of the form @var{f1(xy,t)}.
##@item @var{gD = @{gD1,gD2@}} a cell array with the two components of the prescribed displacements on the boundary section Gamma_1
##@item @var{gN = @{gN1,gN2@}} a cell array with the two components of the surface forces on the boundary section Gamma_2
##@item @var{u0}, @var{v0} the initial displacement and initial velocity.
##@itemize
##@item Any constant function can be given by its scalar value
##@item Any function can be given by a string with the function name
##@item The functions @var{E}, @var{nu} and @var{rho} may also be given as vectors with the values of the function at the Gauss points
##@item The functions @var{f1}, @var{f2}, @var{u0} and @var{v0} may also be given as vectors with the values of the function at the nodes
##@end itemize
##@item @var{t0}, @var{tend} are the initial and final times
##@item @var{steps} is a vector with one or two positive integers.
##@itemize
##@item If @var{steps} = n, then n steps are taken and the n+1 results returned.
##@item If @var{steps} = [n,nint], then n*nint steps are taken and (n+1) results returned.
##@end itemize
##@item @var{options} additional options, given as pairs name/value.
##Currently only the time stepping algorithm can be selected as @var{"solver"} and the possible values
##@itemize
##@item @var{"implicit"} an implicit solver (default)
##@item @var{"explicit"} the standard explicit solver
##@end itemize
##@end itemize
##
##return values
##@itemize
##@item @var{u1}, @var{u2}  matrices with the values of the x- and y-displacements at the nodes. Matrices with n+1 columns with the values of the solution at the nodes at different times @var{t}
##@item @var{t} is the vector with the n+1 values of the times at which the solutions are returned.
##@end itemize
##
## @c Will be cut out in ??? info file and replaced with the same
## @c references explicitly there, since references to core Octave
## @c functions are not automatically transformed from here to there.
## @c BEGIN_CUT_TEXINFO
## @seealso{PlaneStrainDynamic, PlaneStress, PlaneStrain, PlaneStressEig, PlaneStrainEig}
## @c END_CUT_TEXINFO
## @end deftypefn

function [u1,u2,t] = PlaneStressDynamic(Mesh,E,nu,rho,f,gD,gN,u0,v0,t0,tend,steps,varargin)

if ((nargin<12)||(nargin>14))
  print_usage();
endif
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

nGP = size(Mesh.GPT,1); nElem = size(Mesh.elemT,1); nNodes = size(Mesh.nodes,1);

if ischar(f{1})
  f1V = reshape(feval(fFunc{1},Mesh.GP,Mesh.GPT),nGP/nElem,nElem);
elseif isscalar(f{1})
  f1V = f{1}*ones(nGP/nElem,nElem);
else
  f1V = reshape(f{1},nGP/nElem,nElem);
endif

if ischar(f{2})
  f2V = reshape(feval(f{2},Mesh.GP,Mesh.GPT),nGP/nElem,nElem);
elseif isscalar(f{2})
  f2V = f{2}*ones(nGP/nElem,nElem);
else
  f2V = reshape(fFunc{2},nGP/nElem,nElem);
endif

if ischar(u0{1})
  u01V = feval(u0{1},Mesh.nodes);
elseif isscalar(u0{1})
  u01V = u0{1}*ones(nNodes,1);
else
  u01V = u0{1}(:);
endif
if ischar(u0{2})
  u02V = feval(u0{2},Mesh.nodes);
elseif isscalar(u0{2})
  u02V = u0{2}*ones(nNodes,1);
else
  u02V = u0{2}(:);
endif

if ischar(v0{1})
  v01V = feval(v0{1},Mesh.nodes);
elseif isscalar(v0{1})
  v01V = v0{1}*ones(nNodes,1);
else
  v01V = v0{1}(:);
endif
if ischar(v0{2})
  v02V = feval(v0{2},Mesh.nodes);
elseif isscalar(v0{2})
  v02V = v0{2}*ones(nNodes,1);
else
  v02V = v0{2}(:);
endif

%% only solve for the homogeneous BC if necessary
if (gD{1}==0)&&(gD{1}==0)&&(gN{1}==0)&&(gN{1}==0)
  u1_B = zeros(nNodes,1); u2_B = u1_B;
else
  [u1_B,u2_B] = PlaneStress(Mesh,E,nu,{0,0},GD,GN);
endif

switch Mesh.type
  case 'linear'     %% first order elements
    [A,Wrho,Wf] = PStressEquationWM (Mesh,E,nu,rho);
  case 'quadratic'  %% second order elements
    [A,Wrho,Wf] = PStressEquationQuadWM (Mesh,E,nu,rho);
  case 'cubic'      %% third order elements
    [A,Wrho,Wf] = PStressEquationCubicWM (Mesh,E,nu,rho);
endswitch

%%lambda0 = eigs(A,Wrho,1,'sm')  %% smallest eigenvalue
%%lambda0A = eigs(A,1,'sm')      %% smallest eigenvalue
%%lambda0W = eigs(Wrho,1,'sm')   %% smallest eigenvalue
%%PeriodFEM = 2*pi/sqrt(lambda0)

ind_free       = find(Mesh.node2DOF>0);       %% which nodes lead to DOF
ind_Dirichlet  = find(Mesh.node2DOF==0);      %% Dirichlet nodes

t = t0; f1_dep_t = false; f2_dep_t = false;
if ischar(f{1})
  f1Vec = feval(f{1},Mesh.nodes,t);
  f1_dep_t = true;  % has to be evaluated at each timestep
elseif isscalar(f{1})
  f1Vec = f{1}*ones(length(Mesh.nodesT),1);
else
  f1Vec = f{1};
endif

if ischar(f{2})
  f2Vec = feval(f{1},Mesh.nodes,t);
  f2_dep_t = true;  % has to be evaluated at each timestep
elseif isscalar(f{2})
  f2Vec = f{2}*ones(length(Mesh.nodesT),1);
else
  f2Vec = f{2};
endif

if length(steps)==1
  dt = (tend-t0)/steps;
  steps = [steps,1];
else
  dt = (tend-t0)/(steps(1)*steps(2));
endif

switch solver
  case 'IMPLICIT'
    Mleft    = Wrho + dt^2/4*A;   %% matrix for u(t+dt)
    Mmiddle  = 2*Wrho - dt^2/2*A; %% matrix for u(t)
    Mright   = -Mleft;            %% matrix for u(t-dt)
    [L,U,P,Q] = lu(Mleft);        %% P*A*Q = L*U
    t = t0;
    u = zeros(2*nNodes,steps(1)+1);  %% allocate the memory
    u_curr = [u01V-u1_B;u02V-u2_B];  %% current time level
    u_curr = u_curr(ind_free);       %% u(t-dt)
    u_new  = zeros(2*nNodes,1);      %% u(t)
    u(:,1) = [u01V;u02V];
    v0 = [v01V;v02V];
    fVec = [f1Vec;f2Vec]; fVec = fVec(ind_free);

    %% first step
    u_new(ind_free) = Q*(U\(L\(P*(dt*(Wrho+dt^2/4*A)*v0(ind_free)+(Wrho-dt^2/4*A)*u_curr + dt^2/2*Wf*fVec))));
    
    for ii_t = 1:steps(1)
      for ii_2 = 1:steps(2)
	if f1_dep_t
	  f1Vec = feval(f{1},Mesh.nodes,t);
	endif %% f1_dep_t
	if f2_dep_t
	  f2Vec = feval(f{2},Mesh.nodes,t);
	endif %% f2_dep_t
	if or(f1_dep_t,f2_dep_t)
	      fVec = [f1Vec;f2Vec]; fVec = fVec(ind_free);
	endif
	u_temp = Q*(U\(L\(P*(Mmiddle*u_new(ind_free) + Mright*u_curr + dt^2*(Wf*fVec)))));
	u_curr = u_new(ind_free);
	u_new(ind_free) = u_temp;
	t += dt;
      endfor % ii_2
      u(:,ii_t+1) = u_new + [u1_B;u2_B];
    endfor

  case 'EXPLICIT'
    Mleft    = +Wrho;           %% matrix for u(t+dt)
    Mmiddle  = 2*Wrho - dt^2*A; %% matrix for u(t)
    Mright   = -Wrho ;          %% matrix for u(t-dt)
    [L,U,P,Q] = lu(Mleft);  %% P*A*Q = L*U

    lambda = eigs(A,Wrho,1);  %% check for stability
    if(dt>2/sqrt(lambda))
      warning(sprintf(
              'explicit algorithm is unstable, dt = %g, 2/sqrt(lambda) = %g',
              dt,2/sqrt(lambda)))
    endif

    t = t0;
    u = zeros(2*nNodes,steps(1)+1);  %% allocate the memory
    u_curr = [u01V-u1_B;u02V-u2_B];  %% current time level
    u_curr = u_curr(ind_free);       %% u(t-dt)
    u_new  = zeros(2*nNodes,1);      %% u(t)

    u(:,1) = [u01V;u02V];
    v0 = [v01V;v02V];
    fVec = [f1Vec;f2Vec]; fVec = fVec(ind_free);
    u_new(ind_free) = Wrho\( Wrho*dt*v0(ind_free) + Wrho*u_curr
			     + dt^2/2*(Wf*fVec-A*u_curr));    
    for ii_t = 1:steps(1)
      for ii_2 = 1:steps(2)
	if f1_dep_t
	  f1Vec = feval(f{1},Mesh.nodes,t);
	endif %% f1_dep_t
	if f2_dep_t
	  f2Vec = feval(f{2},Mesh.nodes,t);
	endif %% f2_dep_t
	if or(f1_dep_t,f2_dep_t)
	      fVec = [f1Vec;f2Vec]; fVec = fVec(ind_free);
	endif
	u_temp = Q*(U\(L\(P*(Mmiddle*u_new(ind_free) + Mright*u_curr + dt^2*(Wf*fVec)))));
	u_curr = u_new(ind_free);
	u_new(ind_free) = u_temp;
	t += dt;
      endfor % ii_2
      u(:,ii_t+1) = u_new + [u1_B;u2_B];
    endfor

  otherwise
    error('Invalid optional argument for solver: %s, valid are implicit,  explicit',solver);
endswitch

%% dummy output arguments
u1 = u(1:nNodes,:); u2 = u(nNodes+1:end,:); t = linspace(t0,tend,steps(1)+1);
endfunction
