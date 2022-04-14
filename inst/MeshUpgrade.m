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
## @deftypefn{function file}{}@var{MeshQuad} = MeshUpgrade(@var{MeshLin})
##
##   convert a mesh @var{MeshLin} of order 1 to a mesh @var{MeshQuad} of order 2
##
##parameter: @var{MeshLin} the input mesh of order 1
##
##return value: @var{MeshQuad} the output mesh of order 2
##
## @c Will be cut out in ??? info file and replaced with the same
## @c references explicitly there, since references to core Octave
## @c functions are not automatically transformed from here to there.
## @c BEGIN_CUT_TEXINFO
## @seealso{CreateMeshTriangle, CreateMeshRect}
## @c END_CUT_TEXINFO
## @end deftypefn

## Author: Andreas Stahel <andreas.stahel@gmx.com>
## Created: 2020-03-30
## modified for Octave 7.1.0: 2022-04-14

function MeshQuad = MeshUpgrade(Mesh)
  %% MeshQuad = MeshUpgrade(Mesh)
  %% transform a mesh of order 1 to a mesh of order 2

  if nargin ~=1
    print_usage();
  endif
  
  nodes = Mesh.nodes;    nNodes = length(nodes);
  elem  = Mesh.elem;     nElem  = size(elem,1);
  edges = Mesh.edges;
  nodesT= Mesh.nodesT;
  
  ConnMat = sparse(nNodes,nNodes);

  elemQ = zeros(nElem,6);
  for ii = 1:nElem  %% update the elements
    elemQ(ii,1:3) = elem(ii,:);
    %% node 4 in element, between 2 and 3
    newnode = (nodes(elem(ii,2),:) + nodes(elem(ii,3),:))/2;
    ee = [elem(ii,2),elem(ii,3)];
    number =  full(ConnMat(min(ee),max(ee)));
    if number==0  % add a new node
      nNodes = nNodes + 1;
      nodes(nNodes,1:2) = newnode;
      nodesT(nNodes) = Mesh.elemT(ii);
      ConnMat(min(ee),max(ee)) = nNodes;
      elemQ(ii,4) = nNodes;
    else
      elemQ(ii,4) = number;
    endif
    %% node 5 in element, between 1 and 3
    newnode = (nodes(elem(ii,1),:) + nodes(elem(ii,3),:))/2;
    ee = [elem(ii,1),elem(ii,3)];
    number =  full(ConnMat(min(ee),max(ee))); 
    if number==0  % add a new node
      nNodes = nNodes + 1;
      nodes(nNodes,1:2) = newnode;
      nodesT(nNodes) = Mesh.elemT(ii);
      ConnMat(min(ee),max(ee)) = nNodes;
      elemQ(ii,5) = nNodes;
    else
      elemQ(ii,5) = number;
    endif
    %% node 6 in element, between 1 and 2
    newnode = (nodes(elem(ii,1),:) + nodes(elem(ii,2),:))/2;
    ee = [elem(ii,1),elem(ii,2)];
    number =  full(ConnMat(min(ee),max(ee)));
    if number==0  % add a new node
      nNodes = nNodes + 1;
      nodes(nNodes,1:2) = newnode;
      nodesT(nNodes) = Mesh.elemT(ii);
      ConnMat(min(ee),max(ee)) = nNodes;
      elemQ(ii,6) = nNodes;
    else
      elemQ(ii,6) = number;
    endif
  endfor%% update elements
  
  for ii = 1:size(edges,1)  %% update edges
    ee = [edges(ii,1),edges(ii,2)];
    number =  full(ConnMat(min(ee),max(ee)));
%%    nodesT(number) = min([Mesh.nodesT(edges(ii,1)),Mesh.nodesT(edges(ii,2))]);
    nodesT(number) = Mesh.edgesT(ii);
    edges(ii,1:3) =  [edges(ii,1),number,edges(ii,2)];
  endfor %% update edges
  
  MeshQuad.elem   = elemQ;
  MeshQuad.elemT  = Mesh.elemT;
  MeshQuad.elemArea  = Mesh.elemArea;
  MeshQuad.nodes  = nodes;
  MeshQuad.nodesT = nodesT;
  MeshQuad.edges  = edges;
  MeshQuad.edgesT = Mesh.edgesT;

  %% determine area of elements and the GP (Gauss integration Points)
  nElem = size(elemQ)(1);
  GP = zeros(7*nElem,2); MeshQuad.GPT = zeros(7*nElem,1);
  l1 = (12 - 2*sqrt(15))/21;   l2 = (12 + 2*sqrt(15))/21;
  %% coordinates of the Gauss points
  xi = [l1/2, 1-l1, l1/2, l2/2, 1-l2, l2/2, 1/3]';
  nu = [l1/2, l1/2, 1-l1, l2/2, l2/2, 1-l2, 1/3]';

  M = [(1-xi-nu).*(1-2*xi-2*nu) xi.*(2*xi-1) nu.*(2*nu-1) 4*xi.*nu 4*nu.*(1-xi-nu) 4*xi.*(1-xi-nu)];

  %% for each element
  for ne = 1:nElem
    v0 = nodes(elem(ne,1),1:2);
    v1 = nodes(elem(ne,2),1:2)-v0;
    v2 = nodes(elem(ne,3),1:2)-v0;
    GP(7*ne-6:7*ne,:) = v0+xi*v1+nu*v2;
  endfor
  MeshQuad.GP = GP;
  
  ind = (nodesT~=-1);
  MeshQuad.node2DOF = cumsum(ind).*ind;
  MeshQuad.nDOF = sum(ind);
endfunction
