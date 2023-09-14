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

## -*- texinfo -*-
## @deftypefn{function file}{} ShowDeformation(@var{Mesh},@var{u1},@var{u2},@var{factor})
##
##  display the original domain and the deformed domain
##
##parameters:
##@itemize
##@item @var{Mesh} is the mesh describing the domain
##@item @var{u1} vector with the values of the x-displacements at the nodes
##@item @var{u2} vector with the values of the y-displacements at the nodes
##@item @var{factor} scaling factor for the displacements @var{u1} and @var{u2}
##@end itemize
##
## @end deftypefn

## Author: Andreas Stahel <andreas.stahel@gmx.com>
## Created: 2023-09-14

function ShowDeformation(Mesh,u1,u2,factor)
  trimesh(Mesh.elem,Mesh.nodes(:,1),Mesh.nodes(:,2),'color','green','linewidth',1)
   hold on; trimesh(Mesh.elem,Mesh.nodes(:,1)+factor*u1,Mesh.nodes(:,2)+factor*u2,'color','red','linewidth',1)
   hold off
endfunction
