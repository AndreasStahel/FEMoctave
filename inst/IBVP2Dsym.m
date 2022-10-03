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

function [u,t] = IBVP2Dsym(Mesh,m,a,b0,f,gD,gN1,gN2,u0,t0,tend,steps)
  ## -*- texinfo -*-
  ## @deftypefn{function file}{}[@var{u},@var{t}] = IBVP2Dsym(@var{mesh},@var{m},@var{a},@var{b0},@var{f},@var{gD},@var{gN1},@var{gN2},@var{u0},@var{t0},@var{tend},@var{steps})
  ##
  ##   Solve a symmetric initial boundary value problem
  ##
  ##@verbatim
  ## m*d/dt u - div(a*grad u) + b0*u = f         in domain
  ##                               u = gD        on Dirichlet boundary
  ##                    n*(a*grad u) = gN1+gN2*u on Neumann boundary
  ##                           u(t0) = u0        initial value
  ##@end verbatim
  ##
  ##parameters:
  ##@itemize
  ##@item @var{mesh} is the mesh describing the domain and the boundary types
  ##@item @var{m},@var{a},@var{b0},@var{f},@var{gD},@var{gN1},@var{gN2}
  ##are the coefficients and functions describing the PDE.
  ##@*Any constant function can be given by its scalar value.
  ##@*The functions @var{m},@var{a},@var{b0} and @var{f} may also be given as vectors with the values of the function at the Gauss points.
  ##@item @var{f} may be given as a string for a function depending on (x,y) and time t or a a vector with the values at nodes or as scalar.
  ##If @var{f} is given by a scalar or vector it is independent on time.
  ##@item @var{u0} is the initial value, can be given as a constant, function name or as vector with the values at the nodes
  ##@item @var{t0}, @var{tend} are the initial and final times
  ##@item @var{steps} is a vector with one or two positive integers.
  ##@*If @var{steps} = n, then n Crank Nicolson steps are taken and the results returned.
  ##@*If @var{steps} = [n,nint], then n*nint Crank Nicolson steps are taken and (n+1) results returned.
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


if (gD==0)&&(gN1==0)  %% only solve for the homogeneous BC if necessary
  u_B = 0;
else
  u_B = BVP2D(Mesh,a,b0,0,0,0,gD,gN1,0);  %% solve BVP
endif
switch Mesh.type
case 'linear'    %% linear elements
  A = FEMEquation(Mesh,a,b0,0,0,0, 0, 0,gN2);
case 'quadratic' %% quadratic elements
  A = FEMEquationQuad(Mesh,a,b0,0,0,0, 0, 0,gN2);
case 'cubic'     %% cubic elements
  A = FEMEquationCubic(Mesh,a,b0,0,0,0, 0, 0,gN2);
endswitch

Wu  = FEMInterpolWeight(Mesh,m);  %% weight matrix, leading to W* (d/dt u)
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

ind_free = find(Mesh.node2DOF>0);       %% which nodes lead to DOF
ind_Dirichlet = find(Mesh.node2DOF==0); %% Dirichlet nodes
W = Wu(:,ind_free);

t = t0;  f_dep_t = false;
if ischar(f)
  fVec = feval(f,Mesh.nodes,t+dt/2);
  f_dep_t = true;  % has to be evaluated at each timestep
elseif isscalar(f)
  fVec = f*ones(length(Mesh.nodesT),1);
else
  fVec = f;
endif

Mleft = W+dt/2*A;  Mright = W-dt/2*A;
[L,m,Q] = chol(Mleft,'lower');    %% Q'*A*Q = L*L'
Qt = Q'; Lt = L';
t = t0;
u = zeros(length(u0),steps(1)+1);
u_new = u0-u_B;  u_new(ind_Dirichlet) = 0;

u(:,1) = u0;
for ii_t = 1:steps(1)
  for ii_2 = 1:steps(2)
    if f_dep_t
      fVec = feval(u0,Mesh.nodes,t+dt/2);
    endif %% f_dep_t
    u_temp = Q*(Lt\(L\(Qt*(Mright*u_new(ind_free) + dt*(Wf*fVec)))));
    u_new(ind_free) = u_temp;
    t+= dt;
  %    if ischar(gD)  %% evaluate the function
  %      u_new(ind_Dirichlet) = feval(gD,mesh.nodes(ind_Dirichlet,:));
  %    else           %% use the scalar value
  %      u_new(ind_Dirichlet) = gD;
  %    endif %ischar
  endfor % ii_2
  u(:,ii_t+1) = u_new + u_B;
endfor
t = linspace(t0,tend,steps(1)+1);
endfunction
