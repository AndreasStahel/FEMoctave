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
## @deftypefn{function file}{}@var{MeshLin} = MeshQuad2LinearUpgrade(@var{MeshQuad})
##
##   convert a mesh @var{MeshQuad} of order 2 to a mesh @var{MeshLin} of order 1
##
##parameter: @var{MeshQuad} the input mesh of order 2
##
##return value: @var{MeshLin} the output mesh of order 1
##
## @c Will be cut out in ??? info file and replaced with the same
## @c references explicitly there, since references to core Octave
## @c functions are not automatically transformed from here to there.
## @c BEGIN_CUT_TEXINFO
## @seealso{MeshUpgrade, CreateMeshTriangle, CreateMeshRect}
## @c END_CUT_TEXINFO
## @end deftypefn

## Author: Andreas Stahel <andreas.stahel@gmx.com>
## Created: 2020-03-30

function MeshLin = MeshQuad2Linear(MeshQuad)
  %% MeshLin = MeshQuad2Linear(MeshQuad)
  %% transform a mesh of order 1 to a mesh of order 2

  if nargin ~=1
    print_usage();
  endif
  
  nodes = MeshQuad.nodes;    nNodes = length(nodes);
  elem  = MeshQuad.elem;     nElem  = size(elem,1);
  edges = MeshQuad.edges;
  nodesT= MeshQuad.nodesT;
  
  ConnMat = sparse(nNodes,nNodes);

  elemLin = zeros(4*nElem,3);
  elemLinArea = zeros(4*nElem,1);
  elemLinT = zeros(4*nElem,1);
  MeshLin.GPT = zeros(12*nElem,1);
  for ii = 1:nElem  %% update the elements
    elemLin(4*ii-3,:) = elem(ii,[1 6 5]);
    elemLin(4*ii-2,:) = elem(ii,[6 2 4]);
    elemLin(4*ii-1,:) = elem(ii,[4 3 5]);
    elemLin(4*ii  ,:) = elem(ii,[4 5 6]);
    elemLinArea(4*ii-[3 2 1 0]) = MeshQuad.elemArea(ii)/4;
    elemLinT(4*ii-[3 2 1 0])    = MeshQuad.elemT(ii);
    MeshLin.GPT(12*ii -[11:-1:0]) = MeshQuad.GPT(7*ii);
  endfor%% update elements

  edgesLin = zeros(2*size(edges,1),2);
  edgesLinT = zeros(2*size(edges,1),1);
  for ii = 1:size(edges,1)  %% update edges
    edgesLin(2*ii-1,:) = edges(ii,[1,2]);
    edgesLin(2*ii  ,:) = edges(ii,[2 3]);
    edgesLinT(2*ii+[-1,0]) = MeshQuad.edgesT(ii);
  endfor %% update edges
  
  MeshLin.elem   = elemLin;
  MeshLin.elemT  = elemLinT;
  MeshLin.elemArea  =elemLinArea;
  MeshLin.nodes  = nodes;
  MeshLin.nodesT = nodesT;
  MeshLin.edges  = edgesLin;
  MeshLin.edgesT = edgesLinT;
  MeshLin.node2DOF =   MeshQuad.node2DOF;
  MeshLin.nDOF =   MeshQuad.nDOF;

  
  %% determine area of elements and the GP (Gauss integration Points)
  nElem = size(elemLin,1);
  GP = zeros(3*nElem,2); 
  %% for each element
  for ne = 1:nElem
    v0 = nodes(elemLin(ne,1),1:2);
    v1 = nodes(elemLin(ne,2),1:2) - v0;
    v2 = nodes(elemLin(ne,3),1:2) - v0;
    GP(3*ne-2,:) = v0 + v1/6   + v2/6;
    GP(3*ne-1,:) = v0 + v1*2/3 + v2/6;
    GP(3*ne,:)   = v0 + v1/6   + v2*2/3;
    GPT(3*ne-[2 1 0]) = MeshQuad.GPT(ne);
  endfor
  MeshLin.GP = GP;
  MeshLin.GPT = GPT;
endfunction
