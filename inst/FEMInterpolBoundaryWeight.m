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
## @deftypefn{function file}{}@var{WgD} = FEMInterpolBoundaryWeight(@var{FEMmesh},@var{a},@var{b0})
##@*create the matrix to determine the contribution of gD to a IBVP or BVP
##@*the contribution of gD is the determined by WgD*gD, where gD is the vector with the values at the Dirichlet nodes
##
##@verbatim
##     -div(a*grad u)+ b0*u = f         in domain
##                        u = gD        on Dirichlet boundary
##             n*(a*grad u) = gN1+gN2*u on Neumann boundary
##@end verbatim
##
##parameters:
##@itemize
##@item @var{FEMmesh} is the mesh describing the domain and the boundary types.
##@item @var{a},@var{b0} are the coefficients and functions describing the PDE.
##@end itemize
##
##return value:
##@itemize
##@item @var{WgD} is the sparse weight matrix
##@end itemize
##
## @c Will be cut out in ??? info file and replaced with the same
## @c references explicitly there, since references to core Octave
## @c functions are not automatically transformed from here to there.
## @c BEGIN_CUT_TEXINFO
## @seealso{FEMInterpolWeight, IBVP2D, BVP2D, BVP2Dsym, BVP2eig}
## @c END_CUT_TEXINFO
## @end deftypefn

function WgD = FEMInterpolBoundaryWeight(FEMmesh,aFunc,bFunc)
nElem = size(FEMmesh.elem,1); nGP = size(FEMmesh.GP,1);
nEdges = size(FEMmesh.edgesT,1);
DirichletNodes = find (FEMmesh.nodesT==-1);  % numbering of the Dirichlet nodes

if ischar(aFunc)
  aV = reshape(feval(aFunc,FEMmesh.GP,FEMmesh.GPT),nGP/nElem,nElem);
elseif isscalar(aFunc)
  aV = aFunc*ones(nGP/nElem,nElem);
else
  aV = reshape(aFunc,nGP/nElem,nElem);
endif

if ischar(bFunc)
  bV = reshape(feval(bFunc,FEMmesh.GP,FEMmesh.GPT),nGP/nElem,nElem);
elseif isscalar(bFunc)
  bV = bFunc*ones(nGP/nElem,nElem);
else
  bV = reshape(bFunc,nGP/nElem,nElem);
endif

if size(FEMmesh.elem,2)==3  %% linear elements
  %% create memory for the sparse matrix 
  Si = zeros(9*nEdges+9,1); Sj = Si; Sval = Si;% maximal number of contributions
  M = [4 1 1; 1 4 1; 1 1 4]/6;  %% interpolation matrix
  % insert the element matrices and vectors into the global matrix
  ptrDOF = 1;  %% counter for the DOF we are working on
  k1 = [1;2;3]*[1 1 1];
  for k = 1:nElem   %%for each element
    NoNodes = FEMmesh.elem(k,:);
    dofs = FEMmesh.node2DOF(NoNodes);
    if prod(dofs)==0  %% there is a Dirichlet node in this element
      %% compute element stiffness matrix
      cor = FEMmesh.nodes(FEMmesh.elem(k,:),:);  % coordinates of the nodes
      %% compute element stiffness matrix
      area = FEMmesh.elemArea(k);  % area = 0.5*det(T)
      G = [cor(3,2)-cor(2,2),cor(1,2)-cor(3,2),cor(2,2)-cor(1,2);...
	   cor(2,1)-cor(3,1),cor(3,1)-cor(1,1),cor(1,1)-cor(2,1)];
      mat = sum(aV(:,k))/(12*area)*G'*G + area/3*M*diag(bV(:,k))*M;
      for k1 = 1:3
	if dofs(k1)>0
	  for k2 = 1:3
	    if dofs(k2)==0  %% this is a Dirichlet node
	      Si(ptrDOF)   = dofs(k1);    % number of DOF
	      Sj(ptrDOF)   = find(DirichletNodes==NoNodes(k2)); % number of the node
	      Sval(ptrDOF) = mat(k1,k2);
	      ptrDOF++;
	    endif
	  endfor % k2
	endif % dofs(k1)>0
      endfor % k1
    endif % prod(dofs)==0
  endfor % k (elements)
else  %% quadratic elements
  Si = zeros(36*nElem+36,1); Sj=Si; Sval=Si;%% maximal number of contributions

  l1 = (12 - 2*sqrt(15))/21;   l2 = (12 + 2*sqrt(15))/21;
  w1 = (155 - sqrt(15))/2400;  w2 = (155 + sqrt(15))/2400;
  w3 = 0.1125;         w = [w1,w1,w1,w2,w2,w2,w3]';
  xi = [l1/2, 1-l1, l1/2, l2/2, 1-l2, l2/2, 1/3]';
  nu = [l1/2, l1/2, 1-l1, l2/2, l2/2, 1-l2, 1/3]';
  M = [(1-xi-nu).*(1-2*xi-2*nu) xi.*(2*xi-1) nu.*(2*nu-1) 4*xi.*nu 4*nu.*(1-xi-nu) 4*xi.*(1-xi-nu)];
  Mxi = [-3+4*(xi+nu) 4*xi-1 0*xi 4*nu -4*nu  4-8*xi-4*nu];
  Mnu = [-3+4*(xi+nu) 0*xi 4*nu-1 4*xi  4-4*xi-8*nu -4*xi];

  %% insert the element matrices and vectors into the global matrix
  ptrDOF = 1;  %% counter for the DOF we are working on
  k1 = [1:6]'*ones(1,6);
  for k = 1:nElem   %%for each element
    NoNodes = FEMmesh.elem(k,:);
    dofs = FEMmesh.node2DOF(NoNodes);
    if prod(dofs)==0  %% there is a Dirichlet node in this element
      %% compute element stiffness matrix
      cor = FEMmesh.nodes(FEMmesh.elem(k,:),:);  % coordinates of the nodes
      %% compute element stiffness matrix
      area = FEMmesh.elemArea(k);  % area = 0.5*det(T)
      G = [cor(3,2)-cor(2,2),cor(1,2)-cor(3,2),cor(2,2)-cor(1,2);...
	   cor(2,1)-cor(3,1),cor(3,1)-cor(1,1),cor(1,1)-cor(2,1)];
      B = M'*diag(w.*bV(:,k))*M;
      Mtemp = -G(1,2)*Mxi-G(1,3)*Mnu;
      Ax = Mtemp'*diag(w.*aV(:,k))*Mtemp;
      Mtemp = -G(2,2)*Mxi-G(2,3)*Mnu;
      Ay = Mtemp'*diag(w.*aV(:,k))*Mtemp;
      mat = 0.5/area*(Ax+Ay) + 2*area*B;
      for k1 = 1:6
	if dofs(k1)>0
	  for k2 = 1:6
	    if dofs(k2)==0  %% this is a Dirichlet node
	      Si(ptrDOF)   = dofs(k1);    % number of DOF
	      Sj(ptrDOF)   = find(DirichletNodes==NoNodes(k2));%number node
	      Sval(ptrDOF) = mat(k1,k2);
	      ptrDOF++;
	    endif
	  endfor % k2
	endif % dofs(k1)>0
      endfor % k1
    endif % prod(dofs)==0   
  endfor %% k (elements)
endif

Si = Si(1:ptrDOF-1); Sj = Sj(1:ptrDOF-1); Sval = Sval(1:ptrDOF-1);
WgD = sparse(Si,Sj,Sval,FEMmesh.nDOF,length(DirichletNodes));
endfunction
