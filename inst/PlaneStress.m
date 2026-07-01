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


## Author: Andreas Stahel <andreas.stahel@gmx.com>
## Created: 2026-6-25


## -*- texinfo -*-
## @deftypefn{function file}{}[@var{u1},@var{u2}] = PlaneStress(@var{mesh},@var{E},@var{nu},@var{f},@var{gD},@var{gN},@var{options})
##
##   solve a plane stress problem
##
##@verbatim
##    plane stress equation    in domain
##                   u = gD    on Gamma_1
##       force density = gN    on Gamma_2
##       force density = 0     on Gamma_3
##@end verbatim
##
##parameters:
##@itemize
##@item @var{mesh} is the mesh describing the domain and the boundary types
##@item @var{E},@var{nu} Young's modulus and Poisson's ratio for the material
##@item @var{f = @{f1,f2@}} a cell array with the two components of the volume forces
##@item @var{gD = @{gD1,gD2@}} a cell array with the two components of the prescribed displacements on the boundary section Gamma_1
##@item @var{gN = @{gN1,gN2@}} a cell array with the two components of the surface forces on the boundary section Gamma_2
##@item Any constant function can be given by its scalar value
##@item Any function can be given by a string with the function name or a function handle
##@item The functions @var{E}, @var{nu}, @var{f1} and @var{f2} may also be given as vectors with the values of the function at the Gauss points
##@item @var{options} additional options, given as pairs name/value.
##Currently only one option is possible
##@itemize
##@item@var{"thermal"} and value @var{alphaDeltaT} has to be the function to evaluate the product of @var{alpha} (coefficient of heat expansion) and @var{DeltaT} the assumed difference of the temperature. @var{alphaDeltaT} can be given as scalar, vector with the values at the Gauss points or as function or function handle.
##@end itemize
##@end itemize
##
##return values
##@itemize
##@item @var{u1}  vector with the values of the x-displacement at the nodes
##@item @var{u2}  vector with the values of the y-displacement at the nodes
##@end itemize
##
## @c Will be cut out in ??? info file and replaced with the same
## @c references explicitly there, since references to core Octave
## @c functions are not automatically transformed from here to there.
## @c BEGIN_CUT_TEXINFO
## @seealso{PlaneStrain}
## @c END_CUT_TEXINFO
## @end deftypefn

function [u1,u2] = PlaneStress(Mesh,E,nu,f,gD,gN,varargin)
  if nargin < 6
    print_usage();
  endif

  alphaDeltaT = 0;  %% default value
  if (~isempty(varargin))
    for cc = 1:2:length(varargin)
      switch tolower(varargin{cc})
	case {'thermal'}
          alphaDeltaT = varargin{cc+1};
	otherwise
          error('Invalid optional argument, %s',varargin{cc}.name);
      endswitch % switch
    endfor % for
  endif % if

  nElem = size(Mesh.elem,1); nGP = size(Mesh.GP,1);
  function resV = ReadCoefficient(F_Name)
    if ischar(F_Name)
      resV = reshape(feval(F_Name,Mesh.GP,Mesh.GPT),nGP/nElem,nElem);
    elseif is_function_handle(F_Name)
      resV = reshape(F_Name(Mesh.GP),nGP/nElem,nElem);
    elseif isscalar(F_Name)
      resV = F_Name*ones(nGP/nElem,nElem);
    else
      resV = reshape(F_Name,nGP/nElem,nElem);
    endif
  endfunction

  EV            = ReadCoefficient(E);
  nuV           = ReadCoefficient(nu);
  alphaDeltaTV  = ReadCoefficient(alphaDeltaT);
  f1V           = ReadCoefficient(f{1});
  f2V           = ReadCoefficient(f{2});
  ThermalCoeffV = EV.*alphaDeltaTV./(1-nuV);

  switch Mesh.type
    case 'linear'     %% first order elements
      [A,b] = PStressEquationM      (Mesh,EV,nuV,ThermalCoeffV,f1V,f2V,gD,gN);
    case 'quadratic'  %% second order elements
      [A,b] = PStressEquationQuadM  (Mesh,EV,nuV,ThermalCoeffV,f1V,f2V,gD,gN);
    case 'cubic'       %% second order elements
      [A,b] = PStressEquationCubicM (Mesh,EV,nuV,ThermalCoeffV,f1V,f2V,gD,gN);
  endswitch

  ug = -A\b;
  nDOF = Mesh.nDOF;   n = size(Mesh.nodesT,1);  u1 = zeros(n,1); u2 = u1;

  ind_free1 = find(Mesh.node2DOF(:,1)>0);
  ind_free2 = find(Mesh.node2DOF(:,2)>0);

  %%u([ind_free1;n+ind_free2]) = ug;
  u1(ind_free1) = ug([1:nDOF(1)]);
  u2(ind_free2) = ug([nDOF(1)+1:end]);
  ind_Dirichlet1 = find(Mesh.node2DOF(:,1)==0);
  ind_Dirichlet2 = find(Mesh.node2DOF(:,2)==0);

  if is_function_handle(gD{1})
    u1(ind_Dirichlet1)    = gD{1}(Mesh.nodes(ind_Dirichlet1,:));
  elseif ischar(gD{1})
    u1(ind_Dirichlet1)    = feval(gD{1},Mesh.nodes(ind_Dirichlet1,:));
  else u1(ind_Dirichlet1) = gD{1};
  endif

  if is_function_handle(gD{2})
    u2(ind_Dirichlet2)    = gD{2}(Mesh.nodes(ind_Dirichlet2,:));
  elseif ischar(gD{2})
    u2(ind_Dirichlet2)    = feval(gD{2},Mesh.nodes(ind_Dirichlet2,:));
  else u2(ind_Dirichlet2) = gD{2};
  endif
endfunction
