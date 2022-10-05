## Copyright (C) 2022 Andreas Stahel
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn{function file}{}MeshDeformed = MeshDeform(@var{Mesh},@var{Deform})
##
##   Deform the nodes of @var{Mesh} by the transformation @var{Deform}
##
##parameters:
##@itemize
##@item @var{Mesh} the initial mesh with linear elements
##@* this has to be a mesh with linear elements
##@item@var{Deform} the transformation formula
##@* the function @var{Deform} takes one argument @var{xy}, a n by 2 matrix with the x and y components in colums
## and returns the result in a n by 2 matrix.
##@end itemize
##
##return value
##@itemize
##@item @var{DeformedMesh} the deformed mesh consistes of linear elements
##@*use @var{MeshUpgrade} to generate quadratic or cubic elements
##@end itemize
##@end itemize
##
## Sample call:
##@verbatim
## to generate a mesh in polar coordinates:
##function xy_new = Deform(xy)
##  xy_new = [xy(:,1).*cos(xy(:,2)), xy(:,1).*sin(xy(:,2))];
##endfunction
##
##meshDeformed = MeshDeform(mesh,'Deform');
##@end verbatim
## @c Will be cut out in ??? info file and replaced with the same
## @c references explicitly there, since references to core Octave
## @c functions are not automatically transformed from here to there.
## @c BEGIN_CUT_TEXINFO
## @seealso{CreateMeshRect, CreateMeshTriangle, MeshUpgrade, MeshUpgradeCubic}
## @c END_CUT_TEXINFO
## @end deftypefn


## Author: Andreas Stahel <andrea.stahel@gmx.com>
## Created: 2022-09-13

function meshDeformed = MeshDeform (mesh,Deform)

if (nargin ~=2 ) print_usage(); endif

switch mesh.type
  case 'linear'
    meshDeformed = mesh;
    nodes = feval(Deform,mesh.nodes);
    meshDeformed.nodes = nodes;
    for ne = 1:size(mesh.elem,1)
      v0 = nodes(mesh.elem(ne,1),1:2);
      v1 = nodes(mesh.elem(ne,2),1:2) - v0;
      v2 = nodes(mesh.elem(ne,3),1:2) - v0;
      GP(3*ne-2,:) = v0 + v1/6   + v2/6;
      GP(3*ne-1,:) = v0 + v1*2/3 + v2/6;
      GP(3*ne,:)   = v0 + v1/6   + v2*2/3;
      meshDeformed.elemArea(ne) = abs(det([v1;v2]))/2;
    endfor
    meshDeformed.GP = GP;
  otherwise
    disp('MeshDeform() should only be used on meshes with linear elements')
    print_usage()
endswitch

endfunction
