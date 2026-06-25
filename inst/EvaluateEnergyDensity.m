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
##   evaluate the elastic energy density at the nodes for a plane stress or strain setup
##
##parameters:
##@itemize
##@item @var{Mesh} is the mesh describing the domain
##@item @var{eps_xx}, @var{eps_yy}, @var{eps_xy} vectors with the values of the strains at the nodes
##@item @var{E} Young's modulus of elasticity, either as constant or as string with the function name
##@item @var{nu} Poisson's ratio, either as constant or as string with the function name
##@item @var{options} additional options, given as pairs name/value.
##@itemize
##@item @var{"model"} to select the type of plane elasticity problem
##@itemize
##@item @var{"Pstress"} for the plain stress assumption (default)
##@item @var{"Pstrain"} for the plain strain assumption
##@end itemize
##@item @var{"thermal"} and value @var{alphaDeltaT} has to be the function to evaluate the product of @var{alpha} (coefficient of heat expansion) and @var{DeltaT} the assumed difference of the temperature. @var{alphaDeltaT} can be given as scalar, vector with the values at the Gauss points or as function or function handle. Default : @var{alpha} = 0
##@end itemize
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

function w = EvaluateEnergyDensity(Mesh,eps_xx,eps_yy,eps_xy,EFunc,nuFunc,varargin)

Model = 'PSTRESS'; alphaDeltaT = 0;  %% default values
if (~isempty(varargin))
  for cc = 1:2:length(varargin)
    switch tolower(varargin{cc})
      case {'model'}
        Model = toupper(varargin{cc+1});
      case {'thermal'}
        alphaDeltaT = varargin{cc+1};
      otherwise
        error('Invalid optional argument, %s',varargin{cc}.name);
    endswitch % switch
  endfor % for
endif % if


%% evaluate the material parameters at the nodes
function resV = EvaluateFunc(F_name)
  if     ischar(F_nam)              resV = feval(F_name,Mesh.nodes);
  elseif is_function_handle(F_name) resV = F_name(Mesh.nodes);
  else                              resV = F_name*ones(size(Mesh.nodesT,1),1);
  endif
endfunction

EV  = EvaluateFunc(EFunc);
nuV = EvaluateFunc(nuFunc);
alphaDeltaTV = EvaluateFunc(alphaDeltaT);

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
