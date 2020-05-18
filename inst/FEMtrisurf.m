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
## @deftypefn {} {FEMtrisurf (@var{tri}, @var{x}, @var{y}, @var{u})
##
##   display a solution @var{u} as surface on a triangular mesh
##
##parameters:
##@itemize
##@item @var{tri} is the triangulation of the domain
##@*if @var{tri} has three columns a mesh with linear elements is used
##@*if @var{tri} has six columns a mesh with quadratic elements is used
##@item @var{x}, @var{y} coordinates of the nodes in the mesh
##@item @var{u} values of the function to be displayed
##@end itemize
##
## @c Will be cut out in ??? info file and replaced with the same
## @c references explicitly there, since references to core Octave
## @c functions are not automatically transformed from here to there.
## @c BEGIN_CUT_TEXINFO
## @seealso{FEMtrimesh}
## @end deftypefn

## Author: Andreas Stahel <sha1@BiLi>
## Created: 2020-04-15
## the function is a simple wrapper around the standard trisurf

function FEMtrisurf (tri,x,y,u)
  if size(tri,2)==6  % quadratic elements
    tri4 = zeros(4*size(tri,1),3);
    for ii = 1:size(tri)(1)
      tri4(ii*4-3,1:3) = tri(ii,[1 6 5]);
      tri4(ii*4-2,1:3) = tri(ii,[6 2 4]);
      tri4(ii*4-1,1:3) = tri(ii,[5 4 3]);
      tri4(ii*4  ,1:3) = tri(ii,[4 5 6]);
    endfor
  else
    tri4 = tri;
  endif
  trisurf(tri4,x,y,u)
endfunction
