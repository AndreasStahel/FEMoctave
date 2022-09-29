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
## @deftypefn{function file}{}[@var{ux},@var{uy}] = FEMEvaluateGradient(@var{Mesh},@var{u})
##
##   evaluate the gradient of the function u at the nodes
##
##parameters:
##@itemize
##@item @var{Mesh} is the mesh describing the domain and the boundary types
##@item @var{u} vector with the values of the function at the nodes
##@end itemize
##
##return value
##@itemize
##@item @var{ux}  x component of the gradient of u
##@item @var{uy}  y component of the gradient of u
##@end itemize
##@*the values of the gradient are determined on each element
##@*at the nodes the average of the gradients of the neighboring elements is used
##
## @c Will be cut out in ??? info file and replaced with the same
## @c references explicitly there, since references to core Octave
## @c functions are not automatically transformed from here to there.
## @c BEGIN_CUT_TEXINFO
## @seealso{FEMEvalutateGP, BVP2D, BVP2Dsym, BVP2eig, IBVP2D, CreateMeshRect, CreateMeshTriangle}
## @c END_CUT_TEXINFO
## @end deftypefn

## Author: Andreas Stahel <andreas.stahel@gmx.com>
## Created: 2022-09-19

function [ux,uy] = FEMEvaluateGradient(Mesh,u)
% [ux,uy] = FEMEvaluateGradient(Mesh,u)
%  evaluate the gradient [ux,uy] at the nodes
%
% Mesh  is the FEM mesh
% u     vector with the values of the function at the nodes
%
% ux    x component of the gradient of u
% uy    y component of the gradient of u
%
%   the values of the gradient are determined on each element
%   at the nodes the average of the gradients of the elements is used

n_u = length(u);
n_elem = size(Mesh.elemT,1);
ux = zeros(n_u,1); uy = ux; ux_count = ux; % allocate the memory
v1 = Mesh.nodes(Mesh.elem(:,1),:);
v2 = Mesh.nodes(Mesh.elem(:,2),:)-v1;
v3 = Mesh.nodes(Mesh.elem(:,3),:)-v1;

switch Mesh.type
  case 'linear' %% linear elements
    %% slopes in xi and nu direction
    g0   = u(Mesh.elem(:,1));
    g_xi = u(Mesh.elem(:,2)) - g0;
    g_nu = u(Mesh.elem(:,3)) - g0;
    if 0  % loop over elements, eliminated for speed reasons
      for ne = 1:n_elem
	ind = Mesh.elem(ne,:);
	detinv   =  0.5/Mesh.elemArea(ne);
	ux(ind) += ( v3(ne,2)*g_xi(ne) - v2(ne,2)*g_nu(ne))*detinv;
	uy(ind) += (-v3(ne,1)*g_xi(ne) + v2(ne,1)*g_nu(ne))*detinv;
	ux_count(ind) += 1;
      endfor
    else  % no loop solution
      detinv = 0.5./Mesh.elemArea;
      uxs = (+v3(:,2).*g_xi - v2(:,2).*g_nu).*detinv;
      uys = (-v3(:,1).*g_xi + v2(:,1).*g_nu).*detinv;
      for kk = 1:3
	ux(Mesh.elem(:,kk)) += uxs;
	uy(Mesh.elem(:,kk)) += uys;
	ux_count(Mesh.elem(:,kk)) += 1;
      endfor
    endif
    ux = ux./ux_count;
    uy = uy./ux_count;
  case 'quadratic'  %% quadratic elements
    Nxi = [-3 1 1 1 -1 -1;-1 3 -1 1 -1 1; 0 0 0 0 0 0; 0 0 4 2 2 0; 0 0 -4 -2 -2 0;4 -4 0 -2 2 0];
    Nnu = [ -3 1 1 1 -1 -1; 0 0 0 0 0 0; -1 -1 3 1 1 -1;  0 4 0 2 0 2; 4 0 -4 -2 0 2; 0 -4 0 -2 0 -2];
    if 0 %% loop over elements, eliminated for speed reasons
      for ne = 1:n_elem
	ind  = Mesh.elem(ne,:);
	uel = u(ind);  %% values at the nodes
	detinv = 0.5/Mesh.elemArea(ne);
	ux(ind) += ( v3(ne,2)*Nxi'-v2(ne,2)*Nnu')*uel*detinv;
	uy(ind) += (-v3(ne,1)*Nxi'+v2(ne,1)*Nnu')*uel*detinv;
	ux_count(ind) +=1;
      endfor
    else %% no loop
      detinv = 0.5./Mesh.elemArea;
      v2 .*= detinv; v3 .*= detinv;  %% divide by the determinant of T
      uxs = bsxfun(@times,+v3(:,2),u(Mesh.elem)*Nxi)-bsxfun(@times,v2(:,2),u(Mesh.elem)*Nnu);
      uys = bsxfun(@times,-v3(:,1),u(Mesh.elem)*Nxi)+bsxfun(@times,v2(:,1),u(Mesh.elem)*Nnu);
      for kk = 1:6
	ux(Mesh.elem(:,kk)) += uxs(:,kk);
	uy(Mesh.elem(:,kk)) += uys(:,kk);
	ux_count(Mesh.elem(:,kk)) += 1;
      endfor
    endif
    ux = ux./ux_count;
    uy = uy./ux_count;
  case 'cubic'  %% cubic elements
    Nxi = [-11   -2   -2   -2   -2    1   -2   -2    1    1
	           2   11    2    2   -1    2    2   -1    2   -1
        	   0    0    0    0    0    0    0    0    0    0
	           0    0   -9    9    6   -6   -3    0    0    3
	           0    0   18    0    6    6    0    0    0    0
	           0    0  -18    0   -6   -6    0    0    0    0
	           0    0    9    3    6   -6   -9    0    0   -3
	          18    9    0    6    3    0    6   -3   -6   -3
	          -9  -18    0   -6    0   -3   -6    6    3    3
	           0    0    0  -12  -12   12   12    0    0    0]/2;
    Nnu = [-11   -2   -2   -2   -2    1   -2   -2    1    1
	           0    0    0    0    0    0    0    0    0    0
	           2    2   11   -1    2    2   -1    2    2   -1
	           0   18    0    6    0    0    0    0    6    0
        	   0   -9    0    6    9    0    0   -3   -6    3
       	    -9    0  -18    0   -6    3    6   -6   -3    3
	          18    0    9    3    6   -6   -3    6    0   -3
        	   0    9    0    6    3    0    0   -9   -6   -3
        	   0  -18    0   -6    0    0    0    0   -6    0
        	   0    0    0  -12  -12    0    0   12   12    0]/2;
    detinv = 0.5./Mesh.elemArea;
    v2 .*= detinv; v3 .*= detinv;  %% divide by the determinant of T
    uxs = bsxfun(@times,+v3(:,2),u(Mesh.elem)*Nxi)-bsxfun(@times,v2(:,2),u(Mesh.elem)*Nnu);
    uys = bsxfun(@times,-v3(:,1),u(Mesh.elem)*Nxi)+bsxfun(@times,v2(:,1),u(Mesh.elem)*Nnu);
    for kk = 1:10
      ux(Mesh.elem(:,kk)) += uxs(:,kk);
      uy(Mesh.elem(:,kk)) += uys(:,kk);
      ux_count(Mesh.elem(:,kk)) += 1;
    endfor
    ux = ux./ux_count;
    uy = uy./ux_count;
endswitch

endfunction
