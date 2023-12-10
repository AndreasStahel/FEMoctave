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
## @deftypefn{function file}{}@var{mesh} = MeshAddConstraint(@var{Mesh},@var{position},@var{mode})
##
##   apply an additional constraint
##
##parameters:
##@itemize
##@item @var{Mesh} is the mesh describing the domain
##@item @var{position} coordinates of the node, may be approximate
##@item @var{mode} mode of the node with the additional constraint
##@itemize
##@item @var{mode = -1}: fixed value for scalar problems
##@item @var{mode = [-1,-1]}: fixed x and y displacements for elasticity
##@item @var{mode = [-1,-2]}: fixed x-displacement for elasticity
##@item @var{mode = [-2,-1]}: fixed y-displacement for elasticity
##@end itemize
##@end itemize
##
##return values:
##@itemize
##@item @var{mesh} the new mesh with the additional constraints
##@end itemize
##
## @c Will be cut out in ??? info file and replaced with the same
## @c references explicitly there, since references to core Octave
## @c functions are not automatically transformed from here to there.
## @c BEGIN_CUT_TEXINFO
## @seealso{CreateMeshRect, CreateMeshTriangle}
## @c END_CUT_TEXINFO
## @end deftypefn

## Author: Andreas Stahel <andreas.stahel@gmx.com>
## Created: 2023-12-04

function Mesh = MeshAddConstraint(Mesh,Position,Mode)
  %% find the node closest to Position
  [MinVal,MinNo] = min(sum((Mesh.nodes-Position).^2,2));
  Mesh.nodesT(MinNo,:) = Mode; %% apply the Mode
  ind = (Mesh.nodesT~=-1);     %% update the DOF counters
  Mesh.node2DOF = cumsum(ind).*ind;
  Mesh.nDOF = sum(ind);
endfunction
