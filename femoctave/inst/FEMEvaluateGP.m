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
## @deftypefn{function file}{}[@var{uGP},@var{graduGP}] = FEMEvaluateGP(@var{Mesh},@var{u})
##
##   evaluate the function and gradient at the Gauss points
##
##parameters:
##@itemize
##@item @var{Mesh} is the mesh describing the domain and the boundary types
##@item @var{u} vector with the values of at the nodes
##@end itemize
##
##return values
##@itemize
##@item @var{uGP}  values of u at the Gauss points
##@item @var{graduGP}  matrix with the values of the gradient in the columns
##@end itemize
##
## @c Will be cut out in ??? info file and replaced with the same
## @c references explicitly there, since references to core Octave
## @c functions are not automatically transformed from here to there.
## @c BEGIN_CUT_TEXINFO
## @seealso{FEMEvaluateGradient, BVP2D, BVP2Dsym, BVP2eig, IBVP2D, CreateMeshRect, CreateMeshTriangle}
## @c END_CUT_TEXINFO
## @end deftypefn


## Author: Andreas Stahel <andreas.stahel@gmx.com>
## Created: 2020-04-03

function [uGP,graduGP] = FEMEvaluateGP(Mesh,u)
  uGP = zeros(length(Mesh.GP),1);
  graduGP = zeros(length(Mesh.GP),2);

  if size(Mesh.elem,2) == 3  %% linear elements
    Minterp = [4 1 1; 1 4 1; 1 1 4]/6;  %% interpolation matrix
    for ne = 1:length(Mesh.elemT)
      ind = Mesh.elem(ne,:)';
      uGP(ne*3-2:ne*3)= Minterp*u(ind);
      if nargout >=2
	v1   = [Mesh.nodes(ind(1),:),u(ind(1))];
	v2   = [Mesh.nodes(ind(2),:),u(ind(2))]-v1;
	v3   = [Mesh.nodes(ind(3),:),u(ind(3))]-v1;
	grad = cross(v2,v3); grad *= -1/grad(3);
	graduGP(ne*3-2:ne*3,1) = grad(1);
	graduGP(ne*3-2:ne*3,2) = grad(2);
      endif
    endfor%ne
  else  %% quadratic elements
    l1 = (12 - 2*sqrt(15))/21;   l2 = (12 + 2*sqrt(15))/21;
    w1 = (155 - sqrt(15))/2400;  w2 = (155 + sqrt(15))/2400;
    w3 = 0.1125;         w = [w1,w1,w1,w2,w2,w2,w3]';
    
    xi = [l1/2, 1-l1, l1/2, l2/2, 1-l2, l2/2, 1/3]';
    nu = [l1/2, l1/2, 1-l1, l2/2, l2/2, 1-l2, 1/3]';

    %% the interpolation matrices for the function and its partial derivatives
    M = [(1-xi-nu).*(1-2*xi-2*nu) xi.*(2*xi-1) nu.*(2*nu-1) 4*xi.*nu 4*nu.*(1-xi-nu) 4*xi.*(1-xi-nu)];
    Mxi = [-3+4*(xi+nu) 4*xi-1 0*xi 4*nu -4*nu  4-8*xi-4*nu];
    Mnu = [-3+4*(xi+nu) 0*xi 4*nu-1 4*xi  4-4*xi-8*nu -4*xi];
    for ne = 1:length(Mesh.elemT)
      ind = Mesh.elem(ne,:)';
      u_elem = u(ind); % values of u at the nodes
      uGP(ne*7-6:ne*7)= M*u_elem;
      detT = 2*Mesh.elemArea(ne);  % area = 0.5*det(T)
      cor = Mesh.nodes(ind,:);  % coordinates of the nodes
      G = [cor(3,2)-cor(2,2),cor(1,2)-cor(3,2),cor(2,2)-cor(1,2);...
	   cor(2,1)-cor(3,1),cor(3,1)-cor(1,1),cor(1,1)-cor(2,1)];
      graduGP(ne*7-6:ne*7,1) = -(G(1,2)*Mxi+G(1,3)*Mnu)*u_elem/detT;
      graduGP(ne*7-6:ne*7,2) = -(G(2,2)*Mxi+G(2,3)*Mnu)*u_elem/detT;
    endfor %ne
  endif
  
endfunction
