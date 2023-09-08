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
## Created: 2023-09-08


## -*- texinfo -*-
## @deftypefn{function file}{}[@var{la},@var{u1},@var{u2}] = PlaneStressEig(@var{mesh},@var{E},@var{nu},@var{w},@var{nVec},@var{tol})
##
##   solve a plane stress eigenvalue problem
##
##@verbatim
##                 A*u = la*w*u   in domain, plane stress equation
##                   u = 0        on Gamma_1
##       force density = 0        on Gamma_2
##       force density = 0        on Gamma_3
##@end verbatim
##
##parameters:
##@itemize
##@item @var{mesh} is the mesh describing the domain and the boundary types
##@item @var{E},@var{nu} Young's modulus and Poisson's ratio for the material
##@item @var{w} the material density
##@itemize
##@item Any constant function can be given by its scalar value
##@item Any function can be given by a string with the function name
##@item The functions @var{E}, @var{nu} and @var{w}  may also be given as vectors with the values of the function at the Gauss points
##@end itemize
##@item @var{nVec} the number of smallest eigenvalues to be determined
##@item @var{tol} optional tolerance for for the eigenvalue iteration, default 1e-5
##@end itemize
##
##return values
##@itemize
##@item @var{la} the eigenvalues
##@item @var{u1}  matrix with the values of the x-displacement at the nodes
##@item @var{u2}  matrix with the values of the y-displacement at the nodes
##@end itemize
##
## @c Will be cut out in ??? info file and replaced with the same
## @c references explicitly there, since references to core Octave
## @c functions are not automatically transformed from here to there.
## @c BEGIN_CUT_TEXINFO
## @seealso{PlaneStrainEig, PlaneStress, PlaneStrain}
## @c END_CUT_TEXINFO
## @end deftypefn

function [la,u1,u2] = PlaneStressEig(Mesh,E,nu,w,nVec,tol)
  if ((nargin<4)|(nargin>6))
    print_usage();
  endif

  if (nargin==5) tol = 1e-5;endif

  switch Mesh.type
    case 'linear'     %% first order elements
      [A,W] = PStressEquationWM (Mesh,E,nu,w);
    case 'quadratic'  %% second order elements
      [A,W] = PStressEquationQuadWM (Mesh,E,nu,w);
    case 'cubic'      %% third order elements
      [A,W] = PStressEquationCubicWM (Mesh,E,nu,w);
  endswitch

  if nargout == 1
    la = eigSmall(A,W,nVec,tol);
  endif

  if nargout>1
    nDOF = Mesh.nDOF;
    [la,evec] = eigSmall(A,W,nVec,tol);
    ind_free1 = find(Mesh.node2DOF(:,1)>0);
    ind_free2 = find(Mesh.node2DOF(:,2)>0);
    n = size(Mesh.nodes,1);
    u1 = zeros(n,nVec); u2 = u1;
    u1(ind_free1,:) = evec([1:nDOF(1)],:);
    u2(ind_free2,:) = evec([nDOF(1)+1:end],:);
  endif %% nargout >1
endfunction
