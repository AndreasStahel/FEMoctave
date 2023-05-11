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
## @deftypefn{function file}{}[@var{ui},@var{uxi},@var{uyi}] = FEMgriddata(@var{Mesh},@var{u},@var{xi},@var{yi})
##
##   evaluate the function (and gradient) at given points by interpolation
##
##parameters:
##@itemize
##@item @var{Mesh} is the mesh describing the domain
##@* If @var{Mesh} consists of linear elements, piecewise linear interpolation is used.
##@* If @var{Mesh} consists of quadratic elements, piecewise quadratic interpolation is used.
##@* If @var{Mesh} consists of cubic elements, piecewise cubic interpolation is used.
##@item @var{u} vector with the values of the function at the nodes
##@item @var{xi}, @var{yi} coordinates of the points where the function is evaluated
##@end itemize
##
##return values:
##@itemize
##@item @var{ui} values of the interpolated function u
##@item @var{uxi} x component of the gradient of u
##@item @var{uyi} y component of the gradient of u
##@end itemize
##@*The values of the function and the gradient are determined on each element by a piecewise linear, quadratic or cubic interpolation.
##@*If a point is not inside the mesh NaN is returned.
##
## @c Will be cut out in ??? info file and replaced with the same
## @c references explicitly there, since references to core Octave
## @c functions are not automatically transformed from here to there.
## @c BEGIN_CUT_TEXINFO
## @seealso{FEMEvalutateGP, BVP2D, BVP2Dsym, BVP2eig, IBVP2D, CreateMeshRect, CreateMeshTriangle}
## @c END_CUT_TEXINFO
## @end deftypefn


## Author: Andreas Stahel <andreas.stahel@gmx.com>
## Created: 2020-04-03

function [ui,uxi,uyi] = FEMgriddata(Mesh,u,xint,yint)
%  [ui,uxi,uyi] = FEMgriddata(Mesh,u,xint,yint)
%  interpolate the values of the function and its gradient
%
% Mesh  is the FEM mesh
% u     vector with the values of the function at the nodes
% xint  x-coodinates of the points wehre the function is to be evaluated
% yint  y-coodinates of the points wehre the function is to be evaluated
%
% ui    values of the interpolated function
% uxi   x component of the gradient of u
% uyi   y component of the gradient of u

[nr,nc] = size(xint); xint = xint(:); yint = yint(:); NumbPoints = length(xint);
x = Mesh.nodes(:,1);  y = Mesh.nodes(:,2);
ind_elem = tsearch(x,y,Mesh.elem(:,[1:3]),xint,yint);
indNan = find(isnan(ind_elem));  %% points outside of mesh
indOK  = find(~isnan(ind_elem));
ind_elem = ind_elem(indOK);
xint = xint(indOK); yint = yint(indOK);
Mesh.elem = Mesh.elem(ind_elem,:);
v1 = Mesh.nodes(Mesh.elem(:,1),:);
v2 = Mesh.nodes(Mesh.elem(:,2),:)-v1;
v3 = Mesh.nodes(Mesh.elem(:,3),:)-v1;
xint -= v1(:,1); yint -= v1(:,2);
Det = (v2(:,1).*v3(:,2)- v2(:,2).*v3(:,1));
xi  = (v3(:,2).*xint   - v3(:,1).*yint)./Det;
nu  = (v2(:,1).*yint   - v2(:,2).*xint)./Det;

ui = nan(nc*nr,1);       %% default value is NaN
switch Mesh.type
case 'linear'  %% linear elements
  if NumbPoints ==1
    ui(indOK) = [1-xi-nu xi nu]*u(Mesh.elem);
  else
    ui(indOK) = sum([1-xi-nu xi nu].*u(Mesh.elem),2);
  endif
case 'quadratic'  %% quadratic elements
  Onexinu = 1-xi-nu;
  if NumbPoints ==1
    ui(indOK) = [Onexinu.*(1-2*(xi+nu)) xi.*(2*xi-1) nu.*(2*nu-1) 4*xi.*nu 4*nu.*Onexinu 4*xi.*Onexinu]*u(Mesh.elem);
  else
    ui(indOK) = sum([Onexinu.*(1-2*(xi+nu)) xi.*(2*xi-1) nu.*(2*nu-1) 4*xi.*nu 4*nu.*Onexinu 4*xi.*Onexinu].*u(Mesh.elem),2);
  endif
case 'cubic'
  if NumbPoints == 1
  ui(indOK) = [1-11/2*xi-11/2*nu+9*xi.^2+18*xi.*nu+9*nu.^2-9/2*xi.^3-27/2*xi.^2.*nu-27/2*xi.*nu.^2-9/2*nu.^3,...
     xi-9/2*xi.^2+9/2*xi.^3,...
     nu-9/2*nu.^2+9/2*nu.^3,...
     -9/2*xi.*nu+27/2*xi.^2.*nu,...
     -9/2*xi.*nu+27/2*xi.*nu.^2,...
     -9/2*nu+9/2*xi.*nu+18*nu.^2-27/2*xi.*nu.^2-27/2*nu.^3,...
     9*nu-45/2*xi.*nu-45/2*nu.^2+27/2*xi.^2.*nu+27*xi.*nu.^2+27/2*nu.^3,...
     9*xi-45/2*xi.^2-45/2*xi.*nu+27/2*xi.^3+27*xi.^2.*nu+27/2*xi.*nu.^2,...
     -9/2*xi+18*xi.^2+9/2*xi.*nu-27/2*xi.^3-27/2*xi.^2.*nu,...
     27*xi.*nu-27*xi.^2.*nu-27*xi.*nu.^2 ]*u(Mesh.elem);
  else
     ui(indOK) = sum([1-11/2*xi-11/2*nu+9*xi.^2+18*xi.*nu+9*nu.^2-9/2*xi.^3-27/2*xi.^2.*nu-27/2*xi.*nu.^2-9/2*nu.^3,...
     xi-9/2*xi.^2+9/2*xi.^3,...
     nu-9/2*nu.^2+9/2*nu.^3,...
     -9/2*xi.*nu+27/2*xi.^2.*nu,...
     -9/2*xi.*nu+27/2*xi.*nu.^2,...
     -9/2*nu+9/2*xi.*nu+18*nu.^2-27/2*xi.*nu.^2-27/2*nu.^3,...
     9*nu-45/2*xi.*nu-45/2*nu.^2+27/2*xi.^2.*nu+27*xi.*nu.^2+27/2*nu.^3,...
     9*xi-45/2*xi.^2-45/2*xi.*nu+27/2*xi.^3+27*xi.^2.*nu+27/2*xi.*nu.^2,...
     -9/2*xi+18*xi.^2+9/2*xi.*nu-27/2*xi.^3-27/2*xi.^2.*nu,...
     27*xi.*nu-27*xi.^2.*nu-27*xi.*nu.^2 ].*u(Mesh.elem),2);
  endif
endswitch
ui = reshape(ui,nr,nc);  %% return in the desired shape

if nargout > 1  %% determine the gradients too
  uxi = nan(nc*nr,1); uyi = uxi;

  switch Mesh.type
  case 'linear'  %% linear elements
    %% slopes in xi and nu direction
    g_xi = u(Mesh.elem(:,2)) - u(Mesh.elem(:,1));
    g_nu = u(Mesh.elem(:,3)) - u(Mesh.elem(:,1));
  case 'quadratic' %% quadratic elements
    xi4 = 4*xi; nu4 = 4*nu;
    if NumbPoints ==1
      g_xi = [-3+xi4+nu4 xi4-1 zeros(size(xi)) nu4 -nu4 4-2*xi4-nu4]*u(Mesh.elem);
      g_nu = [-3+xi4+nu4 zeros(size(xi)) nu4-1 xi4 4-xi4-2*nu4 -xi4]*u(Mesh.elem);
    else
      g_xi = sum([-3+xi4+nu4 xi4-1 zeros(size(xi)) nu4 -nu4 4-2*xi4-nu4].*u(Mesh.elem),2);
      g_nu = sum([-3+xi4+nu4 zeros(size(xi)) nu4-1 xi4 4-xi4-2*nu4 -xi4].*u(Mesh.elem),2);
    endif
  case 'cubic' %% cubic elements
    z = zeros(size(xi));
    if NumbPoints ==1
      g_xi = [-11/2+18*xi+18*nu-27/2*xi.^2-27*xi.*nu-27/2*nu.^2,...
	      1-9*xi+27/2*xi.^2,...
	      z,...
	      -9/2*nu+27*xi.*nu,...
	      -9/2*nu+27/2*nu.^2,...
	      9/2*nu-27/2*nu.^2,...
	      -45/2*nu+27*xi.*nu+27*nu.^2,...
	      9-45*xi-45/2*nu+81/2*xi.^2+54*xi.*nu+27/2*nu.^2,...
	      -9/2+36*xi+9/2*nu-81/2*xi.^2-27*xi.*nu,...
	      27*nu-54*xi.*nu-27*nu.^2 ]*u(Mesh.elem);
      g_nu = [-11/2+18*xi+18*nu-27/2*xi.^2-27*xi.*nu-27/2*nu.^2,...
	      z,...
	      1-9*nu+27/2*nu.^2,...
	      -9/2*xi+27/2*xi.^2,...
	      -9/2*xi+27*xi.*nu,...
	      -9/2+9/2*xi+36*nu-27*xi.*nu-81/2*nu.^2,...
	      9-45/2*xi-45*nu+27/2*xi.^2+54*xi.*nu+81/2*nu.^2,...
	      -45/2*xi+27*xi.^2+27*xi.*nu,...
	      +9/2*xi-27/2*xi.^2,...
	      27*xi-27*xi.^2-54*xi.*nu ]*u(Mesh.elem);
    else
      g_xi = sum([-11/2+18*xi+18*nu-27/2*xi.^2-27*xi.*nu-27/2*nu.^2,...
		  1-9*xi+27/2*xi.^2,...
		  z,...
		  -9/2*nu+27*xi.*nu,...
		  -9/2*nu+27/2*nu.^2,...
		  9/2*nu-27/2*nu.^2,...
		  -45/2*nu+27*xi.*nu+27*nu.^2,...
		  9-45*xi-45/2*nu+81/2*xi.^2+54*xi.*nu+27/2*nu.^2,...
		  -9/2+36*xi+9/2*nu-81/2*xi.^2-27*xi.*nu,...
		  27*nu-54*xi.*nu-27*nu.^2 ].*u(Mesh.elem),2);
      g_nu = sum([-11/2+18*xi+18*nu-27/2*xi.^2-27*xi.*nu-27/2*nu.^2,...
		  z,...
		  1-9*nu+27/2*nu.^2,...
		  -9/2*xi+27/2*xi.^2,...
		  -9/2*xi+27*xi.*nu,...
		  -9/2+9/2*xi+36*nu-27*xi.*nu-81/2*nu.^2,...
		  9-45/2*xi-45*nu+27/2*xi.^2+54*xi.*nu+81/2*nu.^2,...
		  -45/2*xi+27*xi.^2+27*xi.*nu,...
		  +9/2*xi-27/2*xi.^2,...
		  27*xi-27*xi.^2-54*xi.*nu ].*u(Mesh.elem),2);
    endif
  endswitch
  %%transform to the xy coordinates
  uxi(indOK) = (v3(:,2).*g_xi - v2(:,2).*g_nu)./Det;
  uyi(indOK) = (v2(:,1).*g_nu - v3(:,1).*g_xi)./Det;

  uxi = reshape(uxi,nr,nc);  %% return in the desired shape
  uyi = reshape(uyi,nr,nc);
endif

endfunction
