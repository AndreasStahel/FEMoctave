## Copyright (C) 2022 Andreas Stahel
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
## @deftypefn{function file}{}[@var{sigma_x},@var{sigma_y},@var{tau_xy},@var{sigma_z}] = EvaluateStress(@var{Mesh},@var{u1},@var{u2},@var{E},@var{nu})
##
##   evaluate the normal and shearing stresses at the nodes, using Hooke's law for plane stress or plane strain setups
##
##@itemize
##@item[@var{sigma_x},@var{sigma_y},@var{tau_xy}] = EvaluateStress(@var{Mesh},@var{u1},@var{u2},@var{E},@var{nu})
##@*with three return arguments assumes a plane stress situation
##@item[@var{sigma_x},@var{sigma_y},@var{tau_xy},@var{sigma_z}] = EvaluateStress(@var{Mesh},@var{u1},@var{u2},@var{E},@var{nu})
##@*with four return arguments assumes a plane strain situation
##@end itemize
##
##parameters:
##@itemize
##@item @var{Mesh} is the mesh describing the domain
##@item @var{u1} vector with the values of the x-displacements at the nodes
##@item @var{u2} vector with the values of the y-displacements at the nodes
##@item @var{E} Young's modulus of elasticity, either as constant or as string with the function name
##@item @var{nu} Poisson's ratio, either as constant or as string with the function name
##@end itemize
##
##return values:
##@itemize
##@item @var{sigma_x} values of normal stress in x direction at the nodes
##@item @var{sigma_y} values of normal stress in y direction at the nodes
##@item @var{tau_xy} values of shearing strain at the nodes
##@item @var{sigma_z} values of normal stress in y direction at the nodes, only for plane strain situations
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
## Created: 2022-10-23

function [sigma_x,sigma_y,tau_xy,sigma_z] = EvaluateStress(Mesh,u1,u2,EFunc,nuFunc)
  %% evaluate the strains
  [eps_xx, eps_xy1] = FEMEvaluateGradient(Mesh,u1);
  [eps_xy2,eps_yy]  = FEMEvaluateGradient(Mesh,u2);
  eps_xy = (eps_xy1+eps_xy2)/2;

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

  if nargout == 3 %% plane stress
    %% use Hooke's law
    sigma_x = EV./(1-nuV.^2).*(eps_xx + nuV.*eps_yy);
    sigma_y = EV./(1-nuV.^2).*(nuV.*eps_xx + eps_yy);
    tau_xy  = EV./(1+nuV).*eps_xy;
  endif
  if nargout == 4 %% plane strain
    %% use Hooke's law
    Coeff = EV./((1+nuV).*(1-2*nuV));
    sigma_x = Coeff.*((1-nuV).*eps_xx + nuV.*eps_yy);
    sigma_y = Coeff.*(nuV.*eps_xx + (1-nuV).*eps_yy);
    tau_xy  = Coeff.*(1-2*nuV).*eps_xy;
    sigma_z = Coeff.*nuV.*(eps_xx + eps_yy);
  endif
  
endfunction
