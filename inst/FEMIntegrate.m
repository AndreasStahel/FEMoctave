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
## @deftypefn{function file}{}@var{NumIntegral} = FEMIntegrate(@var{Mesh},@var{u})
##
##   integrate a function u over the domain given in Mesh
##
##parameters:
##@itemize
##@item @var{Mesh} is the mesh describing the domain
##@item @var{u} the function to be integrated
##@*can be given as function name to be evaluated or as scalar value, or as a vector with the values at the nodes or the Gauss points.
##@end itemize
##
##return value
##@itemize
##@item @var{NumIntgeral} the numerical approximation of the integral
##@end itemize
##
## @c Will be cut out in ??? info file and replaced with the same
## @c references explicitly there, since references to core Octave
## @c functions are not automatically transformed from here to there.
## @c BEGIN_CUT_TEXINFO
## @seealso{FEMEvaluateGradient, FEMEvaluateGP, BVP2D, BVP2Dsym, BVP2eig, IBVP2D, CreateMeshRect, CreateMeshTriangle}
## @c END_CUT_TEXINFO
## @end deftypefn


## Author: Andreas Stahel <andreas.stahel@gmx.com>
## Created: 2022-09-19

function NumIntegral = FEMIntegrate(Mesh,u)
  nElem = size(Mesh.elem,1);
  switch Mesh.type
  case 'linear' % % linear elements
    if ischar(u)
      uV = reshape(feval(u,Mesh.GP,Mesh.GPT),3,nElem);
    elseif isscalar(u)
      uV = u*ones(3,nElem);
    elseif length(u) == size(Mesh.nodes,1) %% given at the nodes
      uV = reshape(FEMEvaluateGP(Mesh,u),3,nElem);
    else
      uV = reshape(u,3,nElem);
    endif
    NumIntegral = sum(uV)*(Mesh.elemArea)/3;
  case {'quadratic','cubic'}  %% quadratic or cubic elements
    w1 = (155 - sqrt(15))/2400;  w2 = (155 + sqrt(15))/2400;
    w3 = 0.1125;         w = 2*[w1,w1,w1,w2,w2,w2,w3]';
    if ischar(u)
      uV = reshape(feval(u,Mesh.GP,Mesh.GPT),7,nElem);
    elseif isscalar(u)
      uV = u*ones(7,nElem);
    elseif length(u) == size(Mesh.nodes,1) %% given at the nodes
      uV = reshape(FEMEvaluateGP(Mesh,u),7,nElem);
    else
      uV = reshape(u,7,nElem);
    endif
    NumIntegral = 0;
    for ne = 1:nElem
      NumIntegral += sum(w.*uV(:,ne))*Mesh.elemArea(ne);
    endfor
  endswitch
endfunction
