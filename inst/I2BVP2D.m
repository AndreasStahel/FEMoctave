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

function [u,t] = I2BVP2D(Mesh,m,d,a,b0,bx,by,f,gD,gN1,gN2,u0,v0,t0,tend,steps,varargin)
## -*- texinfo -*-
## @deftypefn{function file}{}[@var{u},@var{t}] = I2BVP2D(@var{Mesh},@var{m},@var{d},@var{a},@var{b0},@var{bx},@var{by},@var{f},@var{gD},@var{gN1},@var{gN2},@var{u0},@var{v0},@var{t0},@var{tend},@var{steps},@var{options})
##
##   Solve an initial boundary value problem
##
##@verbatim
## m*d^2/dt^2 u + 2*d*d/dt u - div(a*grad u-u*(bx,by)) + b0*u = f  in domain
##                        u = gD        on Dirichlet boundary
##  n*(a*grad u -u*(bx,by)) = gN1+gN2*u on Neumann boundary
##                    u(t0) = u0        initial value
##               d/dt u(t0) = v0        initial velocity
##@end verbatim
##
##parameters:
##@itemize
##@item @var{Mesh} is the mesh describing the domain and the boundary types
##@item @var{m},@var{d},@var{a},@var{b0},@var{bx},@var{by},@var{f},@var{gD},@var{gN1},@var{gN2}
##are the coefficients and functions describing the PDE.
##@*Any constant function can be given by its scalar value.
##@*The functions @var{m},@var{d},@var{a},@var{b0},@var{bx},@var{by} and @var{f} may also be given as vectors with the values of the function at the Gauss points.
##@item @var{f} may be given as a string for a function depending on (x,y) and time t or a a vector with the values at nodes or as scalar. If @var{f} is given by a scalar or vector it is independent on time.
##@item @var{u0},@var{v0} are the initial value and velocity, can be given as a constant, function name or as vector with the values at the nodes
##@item @var{t0}, @var{tend} are the initial and final times
##@item @var{steps} is a vector with one or two positive integers.
##@itemize
##@item If @var{steps} = n, then n steps are taken and the n+1 results returned.
##@item If @var{steps} = [n,nint], then n*nint steps are taken and (n+1) results returned.
##@end itemize
##@item @var{options} additional options, given as pairs name/value.
##Currently only the stepping algorithm can be selected as @var{"solver"} and the possible values
##@itemize
##@item @var{"implicit"} the standard implicit solver (default)
##@item @var{"explicit"} the standard explicit solver
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
## @seealso{IBVP2D, BVP2Dsym, BVP2eig, CreateMeshRect, CreateMeshTriangle}
## @c END_CUT_TEXINFO
## @end deftypefn

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

if (gD==0)&&(gN1==0)  %% only solve for the homogeneous BC if necessary
  u_B = 0;
else
  u_B = BVP2D(Mesh,a,b0,bx,by,0,gD,gN1,0);  %% solve BVP
endif
switch Mesh.type
case 'linear'
  A   = FEMEquation(Mesh,a,b0,bx,by,0, 0, 0,gN2);  %% compute with compiled code
case 'quadratic'
  A   = FEMEquationQuad(Mesh,a,b0,bx,by,0, 0, 0,gN2);
case 'cubic'
  A   = FEMEquationCubic(Mesh,a,b0,bx,by,0, 0, 0,gN2);
endswitch

Wf  = FEMInterpolWeight(Mesh,1);  %% weight matrix, leading to W*f
Wu  = FEMInterpolWeight(Mesh,m);  %% weight matrix, leading to W* (d^2/dt^2 u)
D   = FEMInterpolWeight(Mesh,d);  %% weight matrix, leading to D* d/dt u

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

if ischar(v0)
  v0 = feval(v0,Mesh.nodes);
elseif isscalar(u0)
  v0 = v0*ones(length(Mesh.nodesT),1);
else
  v0 = v0(:);
endif

ind_free      = find(Mesh.node2DOF>0);  %% which nodes lead to DOF
ind_Dirichlet = find(Mesh.node2DOF==0); %% Dirichlet nodes
W = Wu(:,ind_free);  D = D (:,ind_free);

t = t0; f_dep_t = false;
if ischar(f)
  fVec = feval(f,Mesh.nodes,t);
  f_dep_t = true;  % has to be evaluated at each timestep
elseif isscalar(f)
  fVec = f*ones(length(Mesh.nodesT),1);
else
  fVec = f;
endif

%%lambda = eigs(A,W,1);
%%disp(sprintf("Values:  lambda = %g, dt = %g, 2/sqrt(lambda) = %g\n",...
%%	      lambda,dt,2/sqrt(lambda)))


switch solver
  case 'IMPLICIT'
    Mleft    = W+dt*D+dt^2/4*A;  %% matrix for u(t+dt)
    Mmiddle  = 2*W - dt^2/2*A;   %% matrix for u(t)
    Mright   = -Mleft + 2*dt*D;  %% matrix for u(t-dt)
    [L,U,P,Q] = lu(Mleft);  %% P*A*Q = L*U
    t = t0;
    u = zeros(length(u0),steps(1)+1);
    
    u_curr = u0-u_B;           %% current time level
    u_curr = u_curr(ind_free); %% u(t-dt)
    u_new  = zeros(size(u0));  %% u(t)
    u(:,1) = u0;

    lambda = eigs(A,W,1);  %% check for stability of first step
    if(dt>2/sqrt(lambda))
      warning(sprintf(
              'first step of implicit algorithm could be unstable, dt = %g, 2/sqrt(lambda) = %g',
              dt,2/sqrt(lambda)))
    endif

    %%    u_new(ind_free) = W\( (W-dt*D)*dt*v0(ind_free) + W*u_curr + dt^2/2*(Wf*fVec-A*u_curr));
    if d==0  %% no damping term
      u_new(ind_free) = Q*(U\(L\(P*( dt*(W+dt^2/4*A)*v0(ind_free) + (W-dt^2/4*A)*u_curr + dt^2/2*Wf*fVec))));
    else
      u_new(ind_free) = (W+dt^2/4*A)\( dt*(W-dt*D+dt^2/4*A)*v0(ind_free) + (W-dt^2/4*A)*u_curr + dt^2/2*Wf*fVec);
    endif
    
    for ii_t = 1:steps(1)
      for ii_2 = 1:steps(2)
	if f_dep_t
	  fVec = feval(f,Mesh.nodes,t);
	endif %% f_dep_t
	u_temp = Q*(U\(L\(P*(Mmiddle*u_new(ind_free) + Mright*u_curr + dt^2*(Wf*fVec)))));
	u_curr = u_new(ind_free);
	u_new(ind_free) = u_temp;
	t+= dt;
      endfor % ii_2
      u(:,ii_t+1) = u_new + u_B;
    endfor
    
  case 'EXPLICIT'
    Mleft    = +W + dt/2*D;  %% matrix for u(t+dt)
    Mmiddle  = 2*W - dt^2*A;   %% matrix for u(t)
    Mright   = -W + dt/2*D;  %% matrix for u(t-dt)
    [L,U,P,Q] = lu(Mleft);  %% P*A*Q = L*U

    lambda = eigs(A,W,1);  %% check for stability
    if(dt>2/sqrt(lambda))
      warning(sprintf(
              'explicit algorithm is unstable, dt = %g, 2/sqrt(lambda) = %g',
              dt,2/sqrt(lambda)))
    endif

    t = t0;
    u = zeros(length(u0),steps(1)+1);
    
    u_curr = u0-u_B;           %% current time level
    u_curr = u_curr(ind_free); %% u(t-dt)
    u_new  = zeros(size(u0));  %% u(t)
    u(:,1) = u0;
    
    u_new(ind_free) = W\( (W-dt*D)*dt*v0(ind_free) + W*u_curr + dt^2/2*(Wf*fVec-A*u_curr));
    for ii_t = 1:steps(1)
      for ii_2 = 1:steps(2)
	if f_dep_t
	  fVec = feval(f,Mesh.nodes,t);
	endif %% f_dep_t
	u_temp = Q*(U\(L\(P*(Mmiddle*u_new(ind_free) + Mright*u_curr + dt^2*(Wf*fVec)))));
	u_curr = u_new(ind_free);
	u_new(ind_free) = u_temp;
	t+= dt;
      endfor % ii_2
      u(:,ii_t+1) = u_new + u_B;
    endfor

  otherwise
    error('Invalid optional argument for solver: %s, valid are implicit,  explicit',solver);
endswitch

t = linspace(t0,tend,steps(1)+1);
endfunction
