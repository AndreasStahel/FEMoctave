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
## @deftypefn{function file}{}@var{FEMmesh} = Delaunay2Mesh(@var{tri},@var{x},@var{y}})
##
##   generate a mesh with elements of order 1, using a Delaunay triangulation
##
##parameters:
##@itemize
##@item @var{tri} the Delaunay triangulation
##@item @var{x},@var{y} the coodinates of the points
##@end itemize
##
##return value
##@itemize
##@item @var{FEMmesh} is the mesh to be used by FEMoctave
##@end itemize
##
## @c END_CUT_TEXINFO
## @end deftypefn

function FEMmesh = Delaunay2Mesh(tri,x,y)
  nElem = size(tri,1);  nNodes = length(x);
  FEMmesh.elemT=ones(nElem,1);
  FEMmesh.nodes=[x(:),y(:)];
  FEMmesh.nodesT=ones(size(x(:)));
  FEMmesh.GP = zeros(3*nElem,2);
  FEMmesh.GPT = -ones(3*nElem,1);
  FEMmesh.elemArea = zeros(nElem,1);
  
  %% detect the edges, thei are only part of one triangle
  ConnMat = sparse(nNodes,nNodes);
  for k = 1:nElem
    v0 = FEMmesh.nodes(tri(k,1),1:2);
    v1 = FEMmesh.nodes(tri(k,2),1:2)-v0;
    v2 = FEMmesh.nodes(tri(k,3),1:2)-v0;
    area = det([v1;v2])/2;
    if area<0  % swap two points to force positive orientation
      tri(k,[2 3]) = tri(k,[3 2]);
      vtemp = v1; v1 = v2; v2 = vtemp;
    endif
    FEMmesh.elemArea(k) = abs(area);
    FEMmesh.GP(3*k-2,:) = v0+v1/6+v2/6;
    FEMmesh.GP(3*k-1,:) = v0+v1*2/3+v2/6;
    FEMmesh.GP(3*k,:)   = v0+v1/6+v2*2/3;
    ind = tri(k,:);
    for ee =[ind(1),ind(2);ind(2),ind(3);ind(3),ind(1)]'
      ConnMat(min(ee),max(ee))++;
    endfor
  endfor
  [ind_i,ind_j] = find( ConnMat==1);
  FEMmesh.edges = [ind_i,ind_j];
  FEMmesh.edgesT = -ones(size(ind_i));
  FEMmesh.nodesT(FEMmesh.edges(:)) = -1;
  ind = (FEMmesh.nodesT ~= -1);
  FEMmesh.node2DOF = cumsum(ind).*ind;
  FEMmesh.nDOF = sum(ind);
  FEMmesh.elem=tri;
endfunction
