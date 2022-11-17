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
## @deftypefn{function file}{}[@var{eps_xx},@var{eps_yy},@var{eps_xy}] = EvaluateStrain(@var{Mesh},@var{u1},@var{u2})
##
##   evaluate the normal and shearing strains at the nodes
##
##parameters:
##@itemize
##@item @var{Mesh} is the mesh describing the domain
##@item @var{u1} vector with the values of the x-displacement at the nodes
##@item @var{u2} vector with the values of the x-displacement's at the nodes
##@end itemize
##
##return values:
##@itemize
##@item @var{eps_xx} values of normal strain in x direction at the nodes
##@item @var{eps_yy} values of normal strain in y direction at the nodes
##@item @var{eps_xy} values of shearing strain at the nodes
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
## Created: 2022-10-23

function [eps_xx,eps_yy,eps_xy] = EvaluateStrain(Mesh,u1,u2)
  [eps_xx, eps_xy1] = FEMEvaluateGradient(Mesh,u1);
  [eps_xy2,eps_yy]  = FEMEvaluateGradient(Mesh,u2);
  eps_xy = (eps_xy1+eps_xy2)/2;
endfunction
