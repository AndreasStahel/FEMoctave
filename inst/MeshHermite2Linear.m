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
## @deftypefn{function file}{}@var{ [MeshLin,uLin]} = MeshHermite2Linear(@var{MeshHermite},@var{uHermite})
##
##   convert a mesh @var{MeshHermite} to a mesh @var{MeshLin} of order 1
##
##parameters:
##@itemize
##@item@var{MeshHermite} the input mesh with Hermite elements
##@item@var{uHermite} the DOF of the solution of the BVP 
##@end itemize
##
##return values:
##@itemize
##@item@var{MeshLin} the output mesh of order 1
##@item@var{uLin} the values of the solution at the nodes
##@end itemize##
## @c Will be cut out in ??? info file and replaced with the same
## @c references explicitly there, since references to core Octave
## @c functions are not automatically transformed from here to there.
## @c BEGIN_CUT_TEXINFO
## @seealso{MeshUpgrade, CreateMeshTriangle, CreateMeshRect, MeshCubic2Linear}
## @c END_CUT_TEXINFO
## @end deftypefn

function [meshL,uL] = MeshHermite2Linear(meshH,uH)
  %% convert a solution uH with Hermite elements to a mesh with linear elements
  elem = meshH.elem; xH = meshH.nodes(:,1); yH = meshH.nodes(:,2);
  n_elem = size(meshH.elem,1);
  n_nodes = (length(meshH.nodesT) - n_elem)/3;
  %% reserve space and save the first entries
  uL = zeros(size(uH)); x = uL; y = uL;
  x(1:n_nodes)  = xH(1:n_nodes);  y(1:n_nodes) = yH(1:n_nodes);
  triL = zeros(3*n_nodes+n_elem,3);
  ConnMat = zeros(n_nodes);  %% is the connection covered?
  counter = n_nodes+1;  %% count the nodes working on
  for jj = 1:n_elem      %% for each element
    NewNodes = zeros(1,7); NewNodes(1:3) = sort(elem(jj,1:3));
    %% add the midpoint
    node = meshH.elem(jj,10); u(counter)  = uH(node);
    x(counter) = xH(node);    y(counter)  = yH(node);
    counter_base = counter;   NewNodes(4) = counter;
    counter = counter+1;
    %% cover the first edge
    ee = NewNodes([1 2]);
    number = ConnMat(ee(1),ee(2));
    if number == 0  %% new edge
      v = [diff(x(ee)),diff(y(ee))];
      x(counter)   = x(ee(1))+  v(1)/3; y(counter)   = y(ee(1))+  v(2)/3;
      x(counter+1) = x(ee(1))+2*v(1)/3; y(counter+1) = y(ee(1))+2*v(2)/3;
      ConnMat(ee(1),ee(2)) = counter;  % mark as worked on
      ConnMat(ee(2),ee(1)) = counter+1;
      NewNodes([5,6]) = counter+[0,1];
      counter += 2;
    else   %% edge already covered
      NewNodes([5,6]) = [number,ConnMat(ee(2),ee(1))];
    endif
    %% cover the second edge
    ee = NewNodes([2 3]);
    number = ConnMat(ee(1),ee(2));
    if number == 0
      v = [diff(x(ee)),diff(y(ee))];
      x(counter)   = x(ee(1))+  v(1)/3; y(counter)   = y(ee(1))+  v(2)/3;
      x(counter+1) = x(ee(1))+2*v(1)/3; y(counter+1) = y(ee(1))+2*v(2)/3;
      ConnMat(ee(1),ee(2)) = counter;  % mark as worked on
      ConnMat(ee(2),ee(1)) = counter+1;
      NewNodes([7 8]) = counter+[0,1];
      counter += 2;
    else   %% edge already covered
      NewNodes([7,8]) = [number,ConnMat(ee(2),ee(1))];
    endif
    %% cover the third edge
    ee = NewNodes([1 3]);
    number = ConnMat(ee(1),ee(2));
    if number == 0
      v = [diff(x(ee)),diff(y(ee))];
      x(counter)   = x(ee(1))+  v(1)/3; y(counter)   = y(ee(1))+  v(2)/3;
      x(counter+1) = x(ee(1))+2*v(1)/3; y(counter+1) = y(ee(1))+2*v(2)/3;
      ConnMat(ee(1),ee(2)) = counter;  % mark as worked on
      ConnMat(ee(2),ee(1)) = counter+1;
      NewNodes([9,10]) = counter+[0,1];
      counter += 2;
    else   %% edge already covered
      NewNodes([9,10]) = [number,ConnMat(ee(2),ee(1))];
    endif
    %% evaluate the points, if required
    if counter>counter_base
      Tri = meshH.nodes(meshH.elem(jj,1:3),:);  % coordinates of the corners
      DOF2BF = fDOF2BF(Tri);
      xt = x(NewNodes)-mean(Tri(:,1));  yt = y(NewNodes)-mean(Tri(:,2));
      BF2uP = [ones(size(xt)),xt,yt,xt.^2,xt.*yt,yt.^2,xt.^3,xt.^2.*yt,xt.*yt.^2,yt.^3];
      uL(NewNodes) =  BF2uP*DOF2BF*uH(meshH.elem(jj,:));
    endif %% counter
    triL((jj-1)*9+[1:9],1:3)= [...
	NewNodes([1,5,9]); NewNodes([5,4,9]);NewNodes([5,6,4]);
	NewNodes([6,7,4]); NewNodes([6,2,7]);NewNodes([9,4,10]);
	NewNodes([4,8,10]);NewNodes([4,7,8]);NewNodes([8,3,10])];
  endfor %%jj for all elements

  nodes = [x,y];
  meshL.elemArea = zeros(9*n_elem,1);
  meshL.GP = zeros(9*n_elem,2);
  for jj = 1:9*n_elem  %% find the Gauss points for the linear elements
    v0 = nodes(triL(jj,1),1:2);
    v1 = nodes(triL(jj,2),1:2)-v0;
    v2 = nodes(triL(jj,3),1:2)-v0;
    area = det([v1;v2])/2;
    if area<0  % swap two points to force positive orientation
      triL(jj,[2 3]) = triL(jj,[3 2]);
      vtemp = v1; v1 = v2; v2 = vtemp;
    endif
    meshL.elemArea(jj) = abs(area);
    meshL.GP(3*jj-2,:) = v0 + v1/6   + v2/6; %% coodinates of the Gauss points
    meshL.GP(3*jj-1,:) = v0 + v1*2/3 + v2/6;
    meshL.GP(3*jj,:)   = v0 + v1/6   + v2*2/3;
  endfor
  meshL.nodes = nodes;
  meshL.elem  = triL;
  meshL.type  = 'linear';
  meshL.elemT = ones(9*n_elem,1);
endfunction
