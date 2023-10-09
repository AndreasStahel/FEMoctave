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
## @deftypefn{function file}{}[@var{eps_xx},@var{eps_yy},@var{eps_zz},@var{eps_xz}] = EvaluateStrainAxi(@var{Mesh},@var{ur},@var{uz})
##
##   evaluate the normal and shearing strains at the nodes
##
##parameters:
##@itemize
##@item @var{Mesh} is the mesh describing the domain
##@item @var{ur} vector with the values of the r-displacements at the nodes
##@item @var{uz} vector with the values of the z-displacements at the nodes
##@end itemize
##
##return values:
##@itemize
##@item @var{eps_xx} values of normal strain in x direction at the nodes
##@item @var{eps_yy} values of normal strain in y direction at the nodes
##@item @var{eps_zz} values of normal strain in z direction at the nodes
##@item @var{eps_xz} values of shearing strain at the nodes
##@end itemize
##
## @c Will be cut out in ??? info file and replaced with the same
## @c references explicitly there, since references to core Octave
## @c functions are not automatically transformed from here to there.
## @c BEGIN_CUT_TEXINFO
## @seealso{EvaluateStress, PlaneStress, PlaneStrain}
## @c END_CUT_TEXINFO
## @end deftypefn

## Author: Andreas Stahel <andreas.stahel@gmx.com>
## Created: 2023-01-08

function [eps_xx,eps_yy,eps_zz,eps_xz] = EvaluateStrainAxi(Mesh,ur,uz)
  [eps_xx, eps_xz1] = FEMEvaluateGradient(Mesh,ur);
  [eps_xz2,eps_zz]  = FEMEvaluateGradient(Mesh,uz);
  eps_xz = (eps_xz1+eps_xz2)/2;
  eps_yy = ur./Mesh.nodes(:,1);    %% has to be fixed at radii 0
  ind = find(Mesh.nodes(:,1)==0);
  eps_yy(ind) = eps_xx(ind);
endfunction
