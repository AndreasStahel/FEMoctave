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
## @deftypefn{function file}{}@var{W} = EvaluateEnergyDensity(@var{Mesh},@var{eps_xx},@var{eps_yy},@var{eps_xy},@var{E},@var{nu},@var{options})
##
##   evaluate the energy density at the nodes for a plane stress or strain setup
##
##parameters:
##@itemize
##@item @var{Mesh} is the mesh describing the domain
##@item @var{eps_xx}, @var{eps_yy}, @var{eps_xy} vectors with the values of the strains at the nodes
##@item @var{E} Young's modulus of elasticity, either as constant or as string with the function name
##@item @var{nu} Poisson's ratio, either as constant or as string with the function name
##@item @var{options} additional options, given as pairs name/value.
##Currently only plain stress or strain be can selected as @var{"model"} and the possible values
##@itemize
##@item @var{"Pstress"} for the plain stress assumption (default)
##@item @var{"Pstrain"} for the plain strain assumption
##@end itemize
##@end itemize
##
##return value:
##@itemize
##@item @var{w} values of the energy density at the nodes
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

function w = EvaluateEnergyDensity(Mesh,eps_xx,eps_yy,eps_xy,EFunc,nuFunc,varargin)

Model = 'PSTRESS';
if (~isempty(varargin))
  for cc = 1:2:length(varargin)
    switch tolower(varargin{cc})
      case {'model'}
        Model = toupper(varargin{cc+1});
      otherwise
        error('Invalid optional argument, %s',varargin{cc}.name);
    endswitch % switch
  endfor % for
endif % if


%% evaluate the material parameters at the nodes
if ischar(EFunc)
  EV = feval(EFunc,Mesh.nodes);
else
  EV = EFunc*ones(size(Mesh.nodesT,1),1);
endif

if ischar(nuFunc)
  nuV = feval(nuFunc,Mesh.nodes);
else
  nuV = nuFunc*ones(size(Mesh.nodesT,1),1);
endif

switch toupper(Model)
case 'PSTRESS'
  w =  EV./(2*(1-nuV.^2)).*(eps_xx.^2+eps_yy.^2+2*nuV.*eps_xx.*eps_yy+2*(1-nuV).*eps_xy.^2);
case 'PSTRAIN'
  EV = EV./(1-nuV.^2);
  nuV = nuV./(1-nuV);
  w =  EV./(2*(1-nuV.^2)).*(eps_xx.^2+eps_yy.^2+2*nuV.*eps_xx.*eps_yy+2*(1-nuV).*eps_xy.^2);
otherwise
  error('Invalid optional argument for model: %s, valid are PSTRESS or PSTRAIN',Model);
endswitch

endfunction
