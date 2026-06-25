## Copyright (C) 2026 Andreas Stahel
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
## @deftypefn{function file}{}[@var{sigma_x},@var{sigma_y},@var{tau_xy},@var{sigma_z}] = EvaluateStress(@var{Mesh},@var{u1},@var{u2},@var{E},@var{nu},@var{options})
##
##   evaluate the normal and shearing stresses at the nodes, using Hooke's law for plane stress or plane strain setups
##
##@itemize
##@item[@var{sigma_x},@var{sigma_y},@var{tau_xy}] = EvaluateStress(@var{Mesh},@var{u1},@var{u2},@var{E},@var{nu},@var{options})
##@*with three return arguments assumes a plane stress situation
##@item[@var{sigma_x},@var{sigma_y},@var{tau_xy},@var{sigma_z}] = EvaluateStress(@var{Mesh},@var{u1},@var{u2},@var{E},@var{nu},@var{options})
##@*with four return arguments assumes a plane strain situation
##@end itemize
##
##parameters:
##@itemize
##@item @var{Mesh} is the mesh describing the domain
##@item @var{u1} vector with the values of the x-displacements at the nodes
##@item @var{u2} vector with the values of the y-displacements at the nodes
##@item @var{E} Young's modulus of elasticity, either as constant or as string with the function name or a function handle
##@item @var{nu} Young's modulus of elasticity, either as constant or as string with the function name or a function handle
##@item @var{options} additional options, given as pairs name/value.
##Currently only one option is possible
##@itemize
##@item@var{"thermal"} and value @var{alphaDeltaT} has to be the function to evaluate the product of @var{alpha} (coefficient of heat expansion) and @var{DeltaT} the assumed difference of the temperature
##@end itemize
##@end itemize
##
##return values:
##@itemize
##@item @var{sigma_x} values of normal stress in x direction at the nodes
##@item @var{sigma_y} values of normal stress in y direction at the nodes
##@item @var{tau_xy} values of shearing strain at the nodes
##@item @var{sigma_z} values of normal stress in z direction at the nodes, only for plane strain situations
##@end itemize
##
## @c Will be cut out in ??? info file and replaced with the same
## @c references explicitly there, since references to core Octave
## @c functions are not automatically transformed from here to there.
## @c BEGIN_CUT_TEXINFO
## @seealso{EvaluateStrain, PlaneStress, PlaneStrain}
## @c END_CUT_TEXINFO
## @end deftypefn

## Author: Andreas Stahel <andreas.stahel@gmx.com>
## Created: 2026-06-22

function [sigma_x,sigma_y,tau_xy,sigma_z] = EvaluateStress(Mesh,u1,u2,EFunc,nuFunc,varargin)

  alphaDeltaT = 0;  %% default value
  if (~isempty(varargin))
    for cc = 1:2:length(varargin)
      switch tolower(varargin{cc})
	case {'thermal'}
          alphaDeltaT = toupper(varargin{cc+1});
	otherwise
          error('Invalid optional argument, %s',varargin{cc}.name);
      endswitch % switch
    endfor % for
  endif % if

  %% evaluate the strains
  [eps_xx, eps_xy1] = FEMEvaluateGradient(Mesh,u1);
  [eps_xy2,eps_yy]  = FEMEvaluateGradient(Mesh,u2);
  eps_xy = (eps_xy1+eps_xy2)/2;

  %% evaluate the material parameters at the nodes
  if ischar(EFunc)                  EV = feval(EFunc,Mesh.nodes);
  elseif is_function_handle(EFunc)  EV = EFunc(Mesh.nodes);
  else                              EV = EFunc*ones(size(Mesh.nodesT,1),1);
  endif

  if ischar(nuFunc)                  nuV = feval(nuFunc,Mesh.nodes);
  elseif  is_function_handle(nuFunc) nuV = feval(nuFunc,Mesh.nodes);
  else                               nuV = nuFunc*ones(size(Mesh.nodesT,1),1);
  endif

  if ischar(alphaDeltaT)
    alphaDeltaTV = feval(alphaDeltaT,Mesh.nodes);
  elseif is_function_handle(alphaDeltaT)
    alphaDeltaTV = alphaDeltaT(Mesh.nodes);
  else
    alphaDeltaTV = alphaDeltaT*ones(size(Mesh.nodesT,1),1);
  endif

  if nargout == 3 %% plane stress
    %% use Hooke's law
    sigma_x = EV./(1-nuV.^2).*(eps_xx+nuV.*eps_yy)-EV.*alphaDeltaTV./(1-nuV);
    sigma_y = EV./(1-nuV.^2).*(nuV.*eps_xx+eps_yy)-EV.*alphaDeltaTV./(1-nuV);
    tau_xy  = EV./(1+nuV).*eps_xy;
  endif
  if nargout == 4 %% plane strain
    %% use Hooke's law
    Coeff = EV./((1+nuV).*(1-2*nuV));
    sigma_x = Coeff.*((1-nuV).*eps_xx + nuV.*eps_yy)-EV.*alphaDeltaTV./(1-2*nuV);
    sigma_y = Coeff.*(nuV.*eps_xx + (1-nuV).*eps_yy)-EV.*alphaDeltaTV./(1-2*nuV);
    tau_xy  = Coeff.*(1-2*nuV).*eps_xy;
    sigma_z = Coeff.*nuV.*(eps_xx + eps_yy)-EV.*alphaDeltaTV./(1-2*nuV);
  endif

  
endfunction
