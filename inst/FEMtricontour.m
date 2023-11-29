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
## @deftypefn {} FEMtricontour (@var{mesh}, @var{u}, @var{v})
##
##   display contours of a function @var{u} on a triangular @var{mesh}
##
##parameters:
##@itemize
##@item @var{mesh} is the mesh
##@item @var{u} values of the function to be displayed
##@item @var{v} contours to be used, default value is 21
##@* if @var{v} is scalar, it is the number of contours
##@* if @var{v} is a vector, it is the levels of the contours
##@end itemize
##
## @c Will be cut out in ??? info file and replaced with the same
## @c references explicitly there, since references to core Octave
## @c functions are not automatically transformed from here to there.
## @c BEGIN_CUT_TEXINFO
## @seealso{FEMtrimesh, FEMtrisurf}
## @end deftypefn

## Author: Andreas Stahel <andreas.stahel@gmx.com>
## Created: 2022-09-13
## the function is a simple wrapper around tricontour

function FEMtricontour (mesh,u,v)
  tri = mesh.elem;
  %% check if there is a contour to draw
  DrawContour = 1;
  if (nargin()==2)||length(v)==1
    if min(u)+eps > max(u)
      DrawContour = 0;
    endif
    else
      if sum((v>min(u)).*(v<max(u)))==0
	DrawContour = 0;
      endif
    endif

  if DrawContour
    switch mesh.type;
      case 'linear'  % linear elements
	triN = tri;
      case 'quadratic'  % quadratic elements
	triN = [tri(:,1),tri(:,6),tri(:,5);
		tri(:,6),tri(:,2),tri(:,4);
		tri(:,5),tri(:,4),tri(:,3);
		tri(:,4),tri(:,5),tri(:,6)];
      case 'cubic'  %% cubic elements
	triN = [tri(:,1), tri(:,8), tri(:,7);
		tri(:,8), tri(:,10),tri(:,7);
		tri(:,8), tri(:,9), tri(:,10);
		tri(:,9), tri(:,4), tri(:,10);
		tri(:,9), tri(:,2), tri(:,4);
		tri(:,7), tri(:,10),tri(:,6);
		tri(:,10),tri(:,5), tri(:,6);
		tri(:,10),tri(:,4), tri(:,5);
		tri(:,6), tri(:,5), tri(:,3)];
    endswitch

    if nargin == 2
      tricontour(triN,mesh.nodes(:,1),mesh.nodes(:,2),u,21);
    else
      tricontour(triN,mesh.nodes(:,1),mesh.nodes(:,2),u,v);
    endif
  else
    warning('FEMtricontour: no contours to be drawn')
  endif %% DrawContour
  
endfunction
