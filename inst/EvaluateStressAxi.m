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
## @deftypefn{function file}{}[@var{sigma_x},@var{sigma_y},@var{sigma_z},@var{tau_xz}] = EvaluateStressAxi(@var{Mesh},@var{ur},@var{uz},@var{E},@var{nu})
##
##   evaluate the normal and shearing stresses at the nodes, using Hooke's law
##
##
##parameters:
##@itemize
##@item @var{Mesh} is the mesh describing the domain
##@item @var{ur} vector with the values of the r-displacements at the nodes
##@item @var{uz} vector with the values of the z-displacements at the nodes
##@item @var{E} Young's modulus of elasticity, either as constant or as string with the function name
##@item @var{nu} Young's modulus of elasticity, either as constant or as string with the function name
##@end itemize
##
##return values:
##@itemize
##@item @var{sigma_x} values of normal stress in x direction at the nodes
##@item @var{sigma_y} values of normal stress in y direction at the nodes
##@item @var{sigma_z} values of normal stress in z direction at the nodes
##@item @var{tau_xz} values of shearing strain at the nodes
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
## Created: 2023-01-08

function [sigma_x,sigma_y,sigma_z,tau_xz] = EvaluateStressAxi(Mesh,ur,uz,EFunc,nuFunc)
  %% evaluate the strains
  [eps_xx, eps_xz1] = FEMEvaluateGradient(Mesh,ur);
  [eps_xz2,eps_zz]  = FEMEvaluateGradient(Mesh,uz);
  eps_yy = ur./Mesh.nodes(:,1);
  eps_xz = (eps_xz1+eps_xz2)/2;

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

  %% use Hooke's law
  Coeff = EV./((1+nuV).*(1-2*nuV));
  sigma_x = Coeff.*((1-nuV).*eps_xx + nuV.*(eps_yy+eps_zz));
  sigma_y = Coeff.*(nuV.*(eps_xx+eps_zz) + (1-nuV).*eps_yy);
  sigma_z = Coeff.*((1-nuV).*eps_zz + nuV.*(eps_xx+eps_yy));
  tau_xz  = Coeff.*(1-2*nuV).*eps_xz;
endfunction
