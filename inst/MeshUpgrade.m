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
## @deftypefn{function file}{}@var{MeshNew} = MeshUpgrade(@var{MeshLin},@var{type})
##
##   convert a mesh @var{MeshLin} of order 1 to a mesh @var{MeshNew} of order 2 or 3
##
##parameters: 
##@itemize
##@item @var{MeshLin} the input mesh of order 1
##@item @var{type} is a string, either 'quadratic' or 'cubic'
##@* the default is 'quadratic'
##@end itemize
##
##return value: @var{MeshNew} the output mesh of order 2 or 3
##
## @c Will be cut out in ??? info file and replaced with the same
## @c references explicitly there, since references to core Octave
## @c functions are not automatically transformed from here to there.
## @c BEGIN_CUT_TEXINFO
## @seealso{CreateMeshTriangle, CreateMeshRect, MeshDeform}
## @c END_CUT_TEXINFO
## @end deftypefn

## Author: Andreas Stahel <andreas.stahel@gmx.com>
## Created: 2022-11-03

function MeshNew = MeshUpgrade(Mesh,type)
  %% MeshNew = MeshUpgrade(Mesh,type)
  %% convert a mesh of order 1 to a mesh of order 2 or 3

  if nargin > 2
    print_usage();
  endif
  if nargin == 1
    type = 'quadratic';
  endif
  
  if ~strcmp(Mesh.type,'linear')
    error("MeshUpgrade: input mesh must consist of linear elements")
  endif
  
  nodes = Mesh.nodes;    nNodes = length(nodes);
  elem  = Mesh.elem;     nElem  = size(elem,1);
  edges = Mesh.edges;
  nodesT= Mesh.nodesT;

  elastic = (Mesh.edgesT(1)<-9);
  ConnMat = zeros(nNodes,nNodes);

  switch tolower(type)
    case "quadratic"
      elemN = zeros(nElem,6);
      for ii = 1:nElem  %% update the elements
	elemN(ii,1:3) = elem(ii,:);
	%% node 4 in element, between 2 and 3
	newnode = (nodes(elem(ii,2),:) + nodes(elem(ii,3),:))/2;
	ee = [elem(ii,2),elem(ii,3)];
	number = ConnMat(min(ee),max(ee));
	if number==0  % add a new node
	  nNodes = nNodes+1;
	  nodes(nNodes,1:2) = newnode;
	  if elastic
	    nodesT(nNodes,[1 2]) = Mesh.elemT(ii)*[1,1];
	  else
	    nodesT(nNodes) = Mesh.elemT(ii);
	  endif
	  ConnMat(min(ee),max(ee)) = nNodes;
	  elemN(ii,4) = nNodes;
	else
	  elemN(ii,4) = number;
	endif
	%% node 5 in element, between 1 and 3
	newnode = (nodes(elem(ii,1),:) + nodes(elem(ii,3),:))/2;
	ee = [elem(ii,1),elem(ii,3)];
	number =  ConnMat(min(ee),max(ee));
	if number==0  % add a new node
	  nNodes = nNodes+1;
	  nodes(nNodes,1:2) = newnode;
	  if elastic
	    nodesT(nNodes,[1 2]) = Mesh.elemT(ii)*[1,1];
	  else
	    nodesT(nNodes) = Mesh.elemT(ii);
	  endif
	  ConnMat(min(ee),max(ee)) = nNodes;
	  elemN(ii,5) = nNodes;
	else
	  elemN(ii,5) = number;
	endif
	%% node 6 in element, between 1 and 2
	newnode = (nodes(elem(ii,1),:) + nodes(elem(ii,2),:))/2;
	ee = [elem(ii,1),elem(ii,2)];
	number =  ConnMat(min(ee),max(ee));
	if number==0  % add a new node
	  nNodes = nNodes+1;
	  nodes(nNodes,1:2) = newnode;
	  if elastic
	    nodesT(nNodes,[1 2]) = Mesh.elemT(ii)*[1,1];
	  else
	    nodesT(nNodes) = Mesh.elemT(ii);
	  endif
	  ConnMat(min(ee),max(ee)) = nNodes;
	  elemN(ii,6) = nNodes;
	else
	  elemN(ii,6) = number;
	endif
      endfor%% update elements

      for ii = 1:size(edges,1)  %% update edges
	ee = [edges(ii,1),edges(ii,2)];
	number =  ConnMat(min(ee),max(ee));
	%%    nodesT(number) = min([Mesh.nodesT(edges(ii,1)),Mesh.nodesT(edges(ii,2))]);
	if elastic
	  nodesT(number,:) = [fix(Mesh.edgesT(ii)/10),mod(Mesh.edgesT(ii),-10)];
	else
	  nodesT(number) = Mesh.edgesT(ii);
	endif
	edges(ii,1:3) =  [edges(ii,1),number,edges(ii,2)];
      endfor %% update edges
      
      MeshNew.type   = 'quadratic';

    case "cubic"
      elemN = zeros(nElem,10);
      for ii = 1:nElem  %% update the elements
	elemN(ii,1:3) = elem(ii,:);
	%% nodes 4 and 5 in element, between 2 and 3
	ee = [elem(ii,2),elem(ii,3)];  
	number = ConnMat(min(ee),max(ee)); % from small index to large index
	direction = (nodes(ee(2),:) - nodes(ee(1),:)); 
	newnode1 = nodes(ee(1),:) + direction/3;
	newnode2 = nodes(ee(1),:) + direction*2/3;
	if number==0  % add two new nodes
	  nNodes = nNodes+1;
	  nodes(nNodes,1:2) = newnode1;
	  if elastic
	    nodesT(nNodes,[1 2]) = Mesh.elemT(ii)*[1,1];
	  else
	    nodesT(nNodes) = Mesh.elemT(ii);
	  endif
	  elemN(ii,4) = nNodes;
	  ConnMat(min(ee),max(ee)) = nNodes;  %% marked the treated edge
	  nNodes = nNodes+1;
	  nodes(nNodes,[1,2])  = newnode2;
	  nodesT(nNodes,[1,2]) = Mesh.elemT(ii);
	  elemN(ii,5) = nNodes;
	else
	  elemN(ii,4) = number+1; elemN(ii,5) = number;  %% swap the order
	endif
	%% nodes 6 and 7 in element, between 3 and 1
	ee = [elem(ii,3),elem(ii,1)];  
	number = ConnMat(min(ee),max(ee)); % from small index to large index
	direction = (nodes(ee(2),:) - nodes(ee(1),:));
	newnode1 = nodes(ee(1),:) + direction/3;
	newnode2 = nodes(ee(1),:) + direction*2/3;
	if number==0  % add two new nodes
	  nNodes = nNodes+1;
	  nodes(nNodes,1:2) = newnode1;
	  if elastic
	    nodesT(nNodes,[1 2]) = Mesh.elemT(ii)*[1,1];
	  else
	    nodesT(nNodes) = Mesh.elemT(ii);
	  endif
	  elemN(ii,6) = nNodes;
	  ConnMat(min(ee),max(ee)) = nNodes;  %% marked the treated edge
	  nNodes = nNodes+1;
	  nodes(nNodes,1:2) = newnode2;
	  nodesT(nNodes,1:2) = Mesh.elemT(ii);
	  elemN(ii,7) = nNodes;
	else
	  elemN(ii,6) = number+1; elemN(ii,7) = number;  %% swap the order
	endif
	%% nodes 8 and 9 in element, between 1 and 2
	ee = [elem(ii,1),elem(ii,2)];  
	number = ConnMat(min(ee),max(ee)); % from small index to large index
	direction = (nodes(ee(2),:) - nodes(ee(1),:)); 
	newnode1 = nodes(ee(1),:) + direction/3;
	newnode2 = nodes(ee(1),:) + direction*2/3;
	if number==0  % add two new nodes
	  nNodes = nNodes+1;
	  nodes(nNodes,1:2) = newnode1;
	  if elastic
	    nodesT(nNodes,[1 2]) = Mesh.elemT(ii)*[1,1];
	  else
	    nodesT(nNodes) = Mesh.elemT(ii);
	  endif
	  nodesT(nNodes) = Mesh.elemT(ii);
	  elemN(ii,8) = nNodes;
	  ConnMat(min(ee),max(ee)) = nNodes;  %% marked the treated edge
	  nNodes = nNodes+1;
	  nodes(nNodes,1:2) = newnode2;
	  nodesT(nNodes,1:2) = Mesh.elemT(ii);
	  elemN(ii,9) = nNodes;
	else
	  elemN(ii,8) = number+1; elemN(ii,9) = number; %% swap the order
	endif
	%% node 10
	nNodes = nNodes+1;
	nodes(nNodes,:) = sum(Mesh.nodes(Mesh.elem(ii,1:3),:))/3;
	nodesT(nNodes,:) = Mesh.elemT(ii);
	elemN(ii,10) = nNodes;
      endfor%% update elements
      
      for ii = 1:size(edges,1)  %% update edges
	ee = [edges(ii,1),edges(ii,2)];
	number =  ConnMat(min(ee),max(ee));
	direction1 = nodes(ee(2),:) - nodes(ee(1),:);
	direction2 = nodes(number+1,:) - nodes(number,:);
	samedir = sign(direction1*direction2');
	%%    nodesT(number) = min([Mesh.nodesT(edges(ii,1)),Mesh.nodesT(edges(ii,2))]);
	if elastic
	  nodesT(number,:)   = [fix(Mesh.edgesT(ii)/10),mod(Mesh.edgesT(ii),-10)];
	  nodesT(number+1,:) = [fix(Mesh.edgesT(ii)/10),mod(Mesh.edgesT(ii),-10)];
	else
	  nodesT(number) = Mesh.edgesT(ii); nodesT(number+1) = Mesh.edgesT(ii);
	endif
	%% one of the two
	if samedir>0
	  edges(ii,1:4) =  [edges(ii,1),number,number+1,edges(ii,2)];
	else
	  edges(ii,1:4) =  [edges(ii,1),number+1,number,edges(ii,2)];
	endif
      endfor %% update edges
      
      MeshNew.type = 'cubic';
    otherwise
      disp(sprintf('the provided type to MeshUpgrade is \"%s\"',type))
      disp(sprintf('permitted are \"quadratic\" or \"cubic\"'))
      print_usage()
  endswitch
  
  MeshNew.elem   = elemN;
  MeshNew.elemT  = Mesh.elemT;
  MeshNew.elemArea  = Mesh.elemArea;
  MeshNew.nodes  = nodes;
  MeshNew.nodesT = nodesT;
  MeshNew.edges  = edges;
  MeshNew.edgesT = Mesh.edgesT;

  %% determine area of elements and the GP (Gauss integration Points)
  nElem = size(elemN)(1);
  GP = zeros(7*nElem,2); MeshNew.GPT = zeros(7*nElem,1);
  l1 = (12 - 2*sqrt(15))/21;   l2 = (12 + 2*sqrt(15))/21;
  %% coordinates of the Gauss points
  xi = [l1/2, 1-l1, l1/2, l2/2, 1-l2, l2/2, 1/3]';
  nu = [l1/2, l1/2, 1-l1, l2/2, l2/2, 1-l2, 1/3]';

  %% for each element
  for ne = 1:nElem
    v0 = nodes(elem(ne,1),1:2);
    v1 = nodes(elem(ne,2),1:2)-v0;
    v2 = nodes(elem(ne,3),1:2)-v0;
    GP(7*ne-6:7*ne,:) = v0+xi*v1+nu*v2;
  endfor
  MeshNew.GP = GP;

  ind = (nodesT~=-1);
  MeshNew.node2DOF = cumsum(ind).*ind;
  MeshNew.nDOF = sum(ind);
endfunction
