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
## @deftypefn{function file}{}@var{wMat} = FEMInterpolWeight(@var{FEMmesh},@var{wFunc})
##@*create the matrix to determine the contribution of w*f to a IBVP or BVP
##@*the contribution of w*f is the determined by wMat*f, where f is the vector with the values at the "free" nodes
##
##@verbatim
##     -div(a*grad u)+ b0*u = w*f       in domain
##                        u = gD        on Dirichlet boundary
##             n*(a*grad u) = gN1+gN2*u on Neumann boundary
##@end verbatim
##
##parameters:
##@itemize
##@item @var{mesh} is the mesh describing the domain and the boundary types
##@item @var{wFunc} is the weight function w
##@*It may be given as a function name, a vector with the values at the Gauss points or as a scalar value
##@end itemize
##
##return value
##@itemize
##@item @var{wMat} is the sparse weight matrix
##@end itemize
##
## @c Will be cut out in ??? info file and replaced with the same
## @c references explicitly there, since references to core Octave
## @c functions are not automatically transformed from here to there.
## @c BEGIN_CUT_TEXINFO
## @seealso{IBVP2D, BVP2D, BVP2Dsym, BVP2eig}
## @c END_CUT_TEXINFO
## @end deftypefn

function [wMat] = FEMInterpolWeight(FEMmesh,wFunc)
nElem = size(FEMmesh.elem,1);  nGP = size(FEMmesh.GP,1);
  
if ischar(wFunc)
  wV = reshape(feval(wFunc,FEMmesh.GP,FEMmesh.GPT),nGP/nElem,nElem);
elseif isscalar(wFunc)
  wV = wFunc*ones(nGP/nElem,nElem);
else
  wV = reshape(wFunc,nGP/nElem,nElem);
endif

if size(FEMmesh.elem,2)==3  %% linear elements
  %% create memory for the sparse matrix 
  Si = zeros(9*nElem+9,1); Sj = Si; Sval = Si;%% maximal number of contributions

  M = [4 1 1; 1 4 1; 1 1 4]/6;  %% interpolation matrix
      % insert the element matrices and vectors into the global matrix
  ptrDOF = 1;  %% counter for the DOF we are working on
  k1 = [1;2;3]*[1 1 1];
  for k = 1:nElem   %%for each element
    %% compute element stiffness matrix and vector
    area = FEMmesh.elemArea(k);
    mat = area/3*M'*diag(wV(:,k))*M;
    NoNodes = FEMmesh.elem(k,:);
    dofs = FEMmesh.node2DOF(NoNodes);
    Si(ptrDOF:ptrDOF+8)   = dofs(k1(:));
    Sj(ptrDOF:ptrDOF+8)   = ([1;1;1]*NoNodes)(:);
    Sval(ptrDOF:ptrDOF+8) = mat(:);
    ptrDOF += 9;
  endfor % k (elements)
else  %% quadratic elements
  Si = zeros(36*nElem+36,1); Sj=Si; Sval=Si;%% maximal number of contributions

  l1 = (12 - 2*sqrt(15))/21;   l2 = (12 + 2*sqrt(15))/21;
  w1 = (155 - sqrt(15))/2400;  w2 = (155 + sqrt(15))/2400;
  w3 = 0.1125;         w = [w1,w1,w1,w2,w2,w2,w3]';
  xi = [l1/2, 1-l1, l1/2, l2/2, 1-l2, l2/2, 1/3]';
  nu = [l1/2, l1/2, 1-l1, l2/2, l2/2, 1-l2, 1/3]';
  M = [(1-xi-nu).*(1-2*xi-2*nu) xi.*(2*xi-1) nu.*(2*nu-1) 4*xi.*nu 4*nu.*(1-xi-nu) 4*xi.*(1-xi-nu)];

  %% insert the element matrices and vectors into the global matrix
  ptrDOF = 1;  %% counter for the DOF we are working on
  k1 = [1:6]'*ones(1,6);
  for k = 1:nElem   %%for each element
    mat  = FEMmesh.elemArea(k)*2*M'*diag(w.*wV(:,k))*M;
    NoNodes = FEMmesh.elem(k,:);
    dofs = FEMmesh.node2DOF(NoNodes);
    Si(ptrDOF:ptrDOF+35)   = dofs(k1(:));
    Sj(ptrDOF:ptrDOF+35)   = (ones(6,1)*NoNodes)(:);
    Sval(ptrDOF:ptrDOF+35) = mat(:);
    ptrDOF += 36;
  endfor %% k (elements)
endif

IND = find((Si.*Sj)>0);
Si = Si(IND); Sj = Sj(IND); Sval = Sval(IND);
wMat = sparse(Si,Sj,Sval,FEMmesh.nDOF,length(FEMmesh.nodesT));
endfunction
