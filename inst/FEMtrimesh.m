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

## -*- texinfo -*- 
## @deftypefn {} FEMtrimesh (@var{mesh}, @var{u})
##
##   display a solution @var{u} on a triangular @var{mesh}
##
##parameters:
##@itemize
##@item @var{mesh} is the mesh
##@item @var{u} values of the function to be displayed
##@* if @var{u} is not given, then the mesh is displayed in 2D
##@end itemize
##
## @c Will be cut out in ??? info file and replaced with the same
## @c references explicitly there, since references to core Octave
## @c functions are not automatically transformed from here to there.
## @c BEGIN_CUT_TEXINFO
## @seealso{FEMtrisurf, FEMtricontour}
## @end deftypefn

## Author: Andreas Stahel <andreas.stahel@gmx.com>
## Created: 2020-04-15
## the function is a simple wrapper around the standard trimesh

function FEMtrimesh (mesh,u)
  x = mesh.nodes(:,1);
  y = mesh.nodes(:,2);
  tri = mesh.elem;
switch mesh.type
  case 'linear'    %% linear elements
    triN = tri;
  case 'quadratic' %% quadratic elements
    triN = [tri(:,1),tri(:,6),tri(:,5);
	    tri(:,6),tri(:,2),tri(:,4);
	    tri(:,5),tri(:,4),tri(:,3);
	    tri(:,4),tri(:,5),tri(:,6)];
  case 'cubic'     %% cubic elements
    triN = [tri(:,1),tri(:,8),tri(:,7);
	    tri(:,8),tri(:,10),tri(:,7);
	    tri(:,8),tri(:,9),tri(:,10);
	    tri(:,9),tri(:,4),tri(:,10);
	    tri(:,9),tri(:,2),tri(:,4);
	    tri(:,7),tri(:,10),tri(:,6);
	    tri(:,10),tri(:,5),tri(:,6);
	    tri(:,10),tri(:,4),tri(:,5);
	    tri(:,6),tri(:,5),tri(:,3)];
  endswitch

  if nargin == 2
    trimesh(triN,mesh.nodes(:,1),mesh.nodes(:,2),u)
  else
    trimesh(triN,mesh.nodes(:,1),mesh.nodes(:,2))
  endif
endfunction
