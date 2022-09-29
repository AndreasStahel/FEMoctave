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

## -*- texinfo -*-
## @deftypefn{function file}{}@var{u} = BVP2Dsym(@var{mesh},@var{a},@var{b0},@var{f},@var{gD},@var{gN1},@var{gN2})
##
##   Solve a symmetric, elliptic boundary value problem
##
##@verbatim
##     -div(a*grad u)+ b0*u = f            in domain
##                        u = gD           on Dirichlet boundary
##             n*(a*grad u) = gN1 + gN2*u  on Neumann boundary
##@end verbatim
##
##parameters:
##@itemize
##@item @var{mesh} is the mesh describing the domain and the boundary types
##@the The element type (linear, qudratic, cubic) is selected by @var{mesh}
##@item @var{a},@var{b0},@var{f},@var{gD},@var{gN1},@var{gN2}
##are the coefficients and functions describing the PDE.
##@*Any constant function can be given by its scalar value.
##@*The functions @var{a},@var{b0} and @var{f} may also be given as vectors
##with the values of the function at the Gauss points.
##@end itemize
##
##return value
##@itemize
##@item @var{u} is the vector with the values of the solution at the nodes
##@end itemize
##
## @c Will be cut out in ??? info file and replaced with the same
## @c references explicitly there, since references to core Octave
## @c functions are not automatically transformed from here to there.
## @c BEGIN_CUT_TEXINFO
## @seealso{BVP2D, BVP2eig, IBVP2D, CreateMeshRect, CreateMeshTriangle}
## @c END_CUT_TEXINFO
## @end deftypefn

%  u = BVP2Dsym(mesh,a,b0,f,gD,gN1,gN2)
%  Solve a symmetric, elliptic boundary value problem
%
%   u = BVP2D(mesh,'a','b0','f','gD','gN1','gN2')
%   u = BVP2D(mesh,aVec,b0Vec,fVec,'gD','gN1','gN2')
%
%      mesh is the mesh describing the domain
%      'a','b0','f','gD','gN1','gN2'
%          are the names of the functions and coefficients
%      the functions 'a','b0' and 'f' may also be given as vector
%          with the values of the function at the Gauss points
%      all function may be given by a scalar value
%
%  -div(a*grad u) + b0*u = f           in domain
%                      u = gD          on Dirichlet boundary
%           n*(a*grad u) = gN1 +gN2*u  on Neumann boundary
%
%
% see also BVP2D, BVP2Deig

function u = BVP2Dsym(Mesh,a,b0,f,gD,gN1,gN2)
  if nargin ~= 7
    print_usage();
  endif

  switch Mesh.type
    case 'linear'  %% first order elements
      [A,b] = FEMEquation (Mesh,a,b0,0,0,f,gD,gN1,gN2); % compute with compiled code
      %% [A,b] = FEMEquationM(Mesh,a,b0,    f,gD,gN1,gN2); % compute with script code
    case 'quadratic' %% second order elements
      %%[A,b] = FEMEquationQuadM(Mesh,a,b0,f,gD,gN1,gN2); % compute with script
      %%A1(find (abs(A)<eps)) = 0; %% remove almost zero elements
      [A,b] = FEMEquationQuad(Mesh,a,b0,0,0,f,gD,gN1,gN2);
  %%A(find (abs(A)<eps)) = 0; %% remove almost zero elements
    case 'cubic' %% third order elements
      %%[As,bs] = FEMEquationCubicM(Mesh,a,b0,f,gD,gN1,gN2);  % compute with script
      [A,b] = FEMEquationCubic(Mesh,a,b0,0,0,f,gD,gN1,gN2); % compute with compiled code
      %%Differences = [ norm(As(:)-A(:)) , norm(bs - b)]
      %%A(find (abs(A)<eps)) = 0; %% remove almost zero elements
  endswitch
  u = FEMSolve(Mesh,A,b,gD);  %% solve the linear system
endfunction
