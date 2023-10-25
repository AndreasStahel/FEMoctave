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


## Author: Andreas Stahel <andreas.stahel@gmx.com>
## Created: 2023-10-26


## -*- texinfo -*-
## @deftypefn{function file}{}[@var{u1},@var{u2},@var{t}] = PlaneStrainDynamic(@var{mesh},@var{E},@var{nu},@var{rho},@var{f},@var{gD},@var{gN},@var{u0},@var{v0},@var{t0},@var{tend},@var{steps},@var{options})
##
##   solve a dynamic plane strain problem
##
##
##parameters:
##@itemize
##@item @var{mesh} is the mesh describing the domain and the boundary types
##@item @var{E},@var{nu} Young's modulus and Poisson's ratio for the material
##@item @var{rho},the density of the material
##@item @var{f = @{f1,f2@}} a cell array with the two components of the volume forces. The functions take two arguments, coordinates xy and time t, i.e. of the form @var{f1(xy,t)}.
##@item @var{gD = @{gD1,gD2@}} a cell array with the two components of the prescribed displacements on the boundary section Gamma_1
##@item @var{gN = @{gN1,gN2@}} a cell array with the two components of the surface forces on the boundary section Gamma_2
##@item @var{u0}, @var{v0} the initial displacement and initial velocity.
##@itemize
##@item Any constant function can be given by its scalar value
##@item Any function can be given by a string with the function name
##@item The functions @var{E}, @var{nu} and @var{rho} may also be given as vectors with the values of the function at the Gauss points
##@item The functions @var{f1}, @var{f2}, @var{u0} and @var{v0} may also be given as vectors with the values of the function at the nodes
##@end itemize
##@item @var{t0}, @var{tend} are the initial and final times
##@item @var{steps} is a vector with one or two positive integers.
##@itemize
##@item If @var{steps} = n, then n steps are taken and the n+1 results returned.
##@item If @var{steps} = [n,nint], then n*nint steps are taken and (n+1) results returned.
##@end itemize
##@item @var{options} additional options, given as pairs name/value.
##Currently only the time stepping algorithm can be selected as @var{"solver"} and the possible values
##@itemize
##@item @var{"implicit"} an implicit solver (default)
##@item @var{"explicit"} the standard explicit solver
##@end itemize
##@end itemize
##
##return values
##@itemize
##@item @var{u1}, @var{u2}  matrices with the values of the x- and y-displacements at the nodes. Matrices with n+1 columns with the values of the solution at the nodes at different times @var{t}
##@item @var{t} is the vector with the n+1 values of the times at which the solutions are returned.
##@end itemize
##
## @c Will be cut out in ??? info file and replaced with the same
## @c references explicitly there, since references to core Octave
## @c functions are not automatically transformed from here to there.
## @c BEGIN_CUT_TEXINFO
## @seealso{PlaneStressDynamic, PlaneStress, PlaneStrain, PlaneStressEig, PlaneStrainEig}
## @c END_CUT_TEXINFO
## @end deftypefn

function [u1,u2,t] = PlaneStrainDynamic(Mesh,E,nu,rho,f,gD,gN,u0,v0,t0,tend,steps,varargin)

if ((nargin<12)||(nargin>14))
  print_usage();
endif
solver = 'IMPLICIT';  %% default value is IMPLICIT

if (~isempty(varargin))
  for cc = 1:2:length(varargin)
    switch tolower(varargin{cc})
      case {'solver'}
        solver = toupper(varargin{cc+1});
      otherwise
        error('Invalid optional argument, %s',varargin{cc}.name);
    endswitch % switch
  endfor % for
endif % if

nElem = size(Mesh.elem,1); nGP = size(Mesh.GP,1);
if ischar(E)
  EV = reshape(feval(E,Mesh.GP,Mesh.GPT),nGP/nElem,nElem);
elseif isscalar(E)
  EV = E*ones(nGP/nElem,nElem);
else
  EV = reshape(E,nGP/nElem,nElem);
endif

if ischar(nu)
  nuV = reshape(feval(nu,Mesh.GP,Mesh.GPT),nGP/nElem,nElem);
elseif isscalar(nu)
  nuV = nu*ones(nGP/nElem,nElem);
else
  nuV = reshape(nu,nGP/nElem,nElem);
endif
Estar = EV./(1-nuV.^2);
nustar = nuV./(1-nuV);

[u1,u2,t] = PlaneStressDynamic(Mesh,Estar,nustar,rho,f,gD,gN,u0,v0,t0,tend,steps,'solver',solver);
endfunction
