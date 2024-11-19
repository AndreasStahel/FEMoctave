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
## @deftypefn{function file}{}@var{u} = BVP2D(@var{mesh},@var{a},@var{b0},@var{bx},@var{by},@var{f},@var{gD},@var{gN1},@var{gN2})
##
##   Solve an elliptic boundary value problem
##
##@verbatim
##  -div(a*grad u - u*(bx,by))+ b0*u = f         in domain
##                                 u = gD        on Dirichlet boundary
##          n*(a*grad u - u*(bx,by)) = gN1+gN2*u on Neumann boundary
##@end verbatim
##
##parameters:
##@itemize
##@item @var{mesh} is the mesh describing the domain and the boundary types
##@item @var{a},@var{b0},@var{bx},@var{by},@var{f},@var{gD},@var{gN1},@var{gN2}
##are the coefficients and functions describing the PDE.
##@*Any constant function can be given by its scalar value.
##@*The functions @var{a},@var{b0},@var{bx},@var{by} and @var{f} may also be given as vectors
##with the values of the function at the Gauss points.
##@item The coefficient @var{a} can also be a symmetric matrix
##A=[axx,axy;axy,ayy] given by the row vector [axx,ayy,axy].
##It can be given as row vector or as string with the function name or
##as nx3 matrix with the values at the Gauss points.

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
## @seealso{BVP2Dsym, BVP2eig, IBVP2D, CreateMeshRect, CreateMeshTriangle}
## @c END_CUT_TEXINFO
## @end deftypefn

function u = BVP2D(Mesh,a,b0,bx,by,f,gD,gN1,gN2)
  if nargin ~=9
    print_usage();
  endif

  switch Mesh.type
    case 'linear' %% first order elements
      [A,b] = FEMEquation    (Mesh,a,b0,bx,by,f,gD,gN1,gN2);
    case 'quadratic'  %% second order elements
      [A,b] = FEMEquationQuad(Mesh,a,b0,bx,by,f,gD,gN1,gN2);
    case 'cubic'  %% third order elements
      [A,b] = FEMEquationCubic(Mesh,a,b0,bx,by,f,gD,gN1,gN2);
  endswitch


  u = FEMSolve(Mesh,A,b,gD);  %% solve the resulting system
endfunction
