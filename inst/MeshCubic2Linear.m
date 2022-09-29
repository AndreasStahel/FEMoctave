## Copyright (C) 2022 Andreas Stahel
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
## @deftypefn{function file}{}@var{MeshLin} = MeshCubic2Linear(@var{MeshCubic})
##
##   convert a mesh @var{MeshCubic} of order 3 to a mesh @var{MeshLin} of order 1
##
##parameter: @var{MeshCubic} the input mesh of order 3
##
##return value: @var{MeshLin} the output mesh of order 1
##
## @c Will be cut out in ??? info file and replaced with the same
## @c references explicitly there, since references to core Octave
## @c functions are not automatically transformed from here to there.
## @c BEGIN_CUT_TEXINFO
## @seealso{MeshUpgrade, CreateMeshTriangle, CreateMeshRect, MeshQuad2Linear}
## @c END_CUT_TEXINFO
## @end deftypefn

## Author: Andreas Stahel <andreas.stahel@gmx.com>
## Created: 2020-03-30

function MeshLin = MeshCubic2Linear(MeshCubic)
  %% MeshLin = MeshCubic2Linear(MeshCubic)
  %% transform a mesh of order 3 to a mesh of order 1

  if nargin ~=1
    print_usage();
  endif
  
  nodes = MeshCubic.nodes;    nNodes = length(nodes);
  elem  = MeshCubic.elem;     nElem  = size(elem,1);
  edges = MeshCubic.edges;
  nodesT= MeshCubic.nodesT;
  
  ConnMat = sparse(nNodes,nNodes);

  elemLin = zeros(9*nElem,3);
  elemLinArea = zeros(9*nElem,1);
  elemLinT = zeros(9*nElem,1);
  MeshLin.GPT = zeros(9*3*nElem,1);
  for ii = 1:nElem  %% update the elements
    elemLin(9*ii-8,:) = elem(ii,[1 8 7]);
    elemLin(9*ii-7,:) = elem(ii,[8 10 7]);
    elemLin(9*ii-6,:) = elem(ii,[8 9 10]);
    elemLin(9*ii-5,:) = elem(ii,[9 4 10]);
    elemLin(9*ii-4,:) = elem(ii,[9 2 4]);
    elemLin(9*ii-3,:) = elem(ii,[7 10 6]);
    elemLin(9*ii-2,:) = elem(ii,[10 5 6]);
    elemLin(9*ii-1,:) = elem(ii,[10 4 5]);
    elemLin(9*ii  ,:) = elem(ii,[6 5 3]);
    elemLinArea(9*ii-[8:-1:0]) = MeshCubic.elemArea(ii)/9;
    elemLinT(9*ii-[8:-1:0])    = MeshCubic.elemT(ii);
    MeshLin.GPT(9*3*ii -[9*3-1:-1:0]) = MeshCubic.GPT(7*ii);
  endfor%% update elements

  edgesLin  = zeros(3*size(edges,1),2);
  edgesLinT = zeros(3*size(edges,1),1);
  for ii = 1:size(edges,1)  %% update edges
    edgesLin(3*ii-2,:) = edges(ii,[1,2]);
    edgesLin(3*ii-1,:) = edges(ii,[2,3]);
    edgesLin(3*ii  ,:) = edges(ii,[3 4]);
    edgesLinT(3*ii+[-2,-1,0]) = MeshCubic.edgesT(ii);
  endfor %% update edges
  
  MeshLin.elem      = elemLin;
  MeshLin.elemT     = elemLinT;
  MeshLin.elemArea  = elemLinArea;
  MeshLin.nodes     = nodes;
  MeshLin.nodesT    = nodesT;
  MeshLin.edges     = edgesLin;
  MeshLin.edgesT    = edgesLinT;
  MeshLin.node2DOF  =   MeshCubic.node2DOF;
  MeshLin.nDOF      =   MeshCubic.nDOF;
  MeshLin.type      = 'linear';
  
  %% determine area of elements and the GP (Gauss integration Points)
  nElem = size(elemLin,1);
  GP    = zeros(3*nElem,2); 
  %% for each element
  for ne = 1:nElem
    v0 = nodes(elemLin(ne,1),1:2);
    v1 = nodes(elemLin(ne,2),1:2) - v0;
    v2 = nodes(elemLin(ne,3),1:2) - v0;
    GP(3*ne-2,:) = v0 + v1/6   + v2/6;
    GP(3*ne-1,:) = v0 + v1*2/3 + v2/6;
    GP(3*ne,:)   = v0 + v1/6   + v2*2/3;
  endfor
  MeshLin.GP  = GP;
endfunction
