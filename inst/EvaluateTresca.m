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
## @deftypefn{function file}{}@var{Tresca} = EvaluateTresca(@var{sigma_x},@var{sigma_y},@var{tau_xy},@var{sigma_z})
##
##   evaluate the Tresca stress at the nodes
##
##@itemize
##@item @var{Tresca} = EvaluateTresca(@var{sigma_x},@var{sigma_y},@var{tau_xy})
##@*with three input arguments assumes a plane stress situation
##@item @var{Tresca} = EvaluateTresca(@var{sigma_x},@var{sigma_y},@var{tau_xy},@var{sigma_z})
##@*with four input arguments assumes a plane strain situation
##@end itemize
##
##parameters:
##@itemize
##@item @var{sigma_x} values of normal stress in x direction at the nodes
##@item @var{sigma_y} values of normal stress in y direction at the nodes
##@item @var{tau_xy} values of shearing strain at the nodes
##@item @var{sigma_z} values of normal stress in z direction at the nodes, only for plane strain situations
##@end itemize
##
##return values:
##@itemize
##@item @var{Tresca} Tresca stress at the nodes
##@end itemize
##
## @c Will be cut out in ??? info file and replaced with the same
## @c references explicitly there, since references to core Octave
## @c functions are not automatically transformed from here to there.
## @c BEGIN_CUT_TEXINFO
## @seealso{EvaluateStress, EvaluateStrain, PlaneStress, PlaneStrain}
## @c END_CUT_TEXINFO
## @end deftypefn

## Author: Andreas Stahel <andreas.stahel@gmx.com>
## Created: 2022-10-23

function sigma_T = EvaluateTresca(sigma_x,sigma_y,tau_xy,varargin)
  Discr = sqrt((sigma_x-sigma_y).^2+4*tau_xy.^2)/2;
  Sum = sigma_x+sigma_y;
  sigma_1 = Sum + Discr; sigma_2 = Sum - Discr;
  if nargin == 3  %% plane stress
    sigma_T = max([abs(2*Discr),abs(sigma_1),abs(sigma_2)],[],2);
  else  %% plane stress
    sigma_3 = varargin{1};
    sigma_T = max([abs(sigma_1-sigma_2),abs(sigma_1-sigma_3),abs(sigma_2-sigma_3)],[],2);
  endif
  
  endfunction

