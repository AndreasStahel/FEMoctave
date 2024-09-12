## Copyright (C) 2023 Andreas Stahel
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
## @deftypefn{function file}{}@var{W} = EvaluateEnergyDensityAxi(@var{Mesh},@var{eps_xx},@var{eps_yy},@var{eps_zz},@var{eps_xz},@var{E},@var{nu})
##
##   evaluate the elastic energy density at the nodes for an axially symmetric setup
##
##parameters:
##@itemize
##@item @var{Mesh} is the mesh describing the domain
##@item @var{eps_xx}, @var{eps_yy}, @var{eps_zz}, @var{eps_xz} vectors with the values of the strains at the nodes in the plane y=0
##@item @var{E} Young's modulus of elasticity, either as constant or as string with the function name
##@item @var{nu} Poisson's ratio, either as constant or as string with the function name
##@end itemize
##
##return value:
##@itemize
##@item @var{w} values of the elastic energy density at the nodes
##@end itemize
##
## @c Will be cut out in ??? info file and replaced with the same
## @c references explicitly there, since references to core Octave
## @c functions are not automatically transformed from here to there.
## @c BEGIN_CUT_TEXINFO
## @seealso{EvaluateStrain, EvaluateStress, EvaluateVonMises, EvaluateTresca, EvaluatePrincipalStress}
## @c END_CUT_TEXINFO
## @end deftypefn

## Author: Andreas Stahel <andreas.stahel@gmx.com>
## Created: 2022-11-04

function w = EvaluateEnergyDensityAxi(Mesh,eps_xx,eps_yy,eps_zz,eps_xz,EFunc,nuFunc)

%% evaluate the material parameters at the nodes
if ischar(EFunc)
  EV = feval(EFunc,Mesh.nodes);
else
  %%EV = EFunc*ones(size(Mesh.nodesT,1),1);
  EV = EFunc;
endif

if ischar(nuFunc)
  nuV = feval(nuFunc,Mesh.nodes);
else
  %%nuV = nuFunc*ones(size(Mesh.nodesT,1),1);
  nuV = nuFunc;
endif

w =  EV./(2*(1-2*nuV).*(1+nuV)).*((1-nuV).*(eps_xx.^2+eps_yy.^2+eps_zz.^2)+...
      2*nuV.*(eps_xx.*eps_yy + eps_yy.*eps_zz + eps_zz.*eps_xx))...
     + EV./(1+nuV).*(eps_xz.^2);
endfunction
