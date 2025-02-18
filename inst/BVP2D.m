## Copyright (C) 2025 Andreas Stahel
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
## Created: 2025-01-04


## -*- texinfo -*-
## @deftypefn{function file }{}@var{ u} = BVP2D(@var{mesh},@var{a},@var{b0},@var{bx},@var{by},@var{f},@var{gD},@var{gN1},@var{gN2},@var{options})
##
##   Solve an elliptic boundary value problem
##
##@verbatim
##  -div(a*grad u - u*(bx,by))+ b0*u = f         in domain
##                                 u = gD        on Dirichlet boundary
##          n*(a*grad u - u*(bx,by)) = gN1+gN2*u on Neumann boundary
##@end verbatim
##
##parameters:
##@itemize
##@item @var{mesh} is the mesh describing the domain and the boundary types
##@item @var{a},@var{b0},@var{bx},@var{by},@var{f},@var{gD},@var{gN1},@var{gN2}
##are the coefficients and functions describing the PDE.
##@*Any constant function can be given by its scalar value.
##@*The functions @var{a},@var{b0},@var{bx},@var{by} and @var{f} may also be given as vectors
##with the values of the function at the Gauss points.
##@item The coefficient @var{a} can also be a symmetric matrix
##A=[axx,axy;axy,ayy] given by the row vector [axx,ayy,axy].
##It can be given as row vector or as string with the function name or
##as nx3 matrix with the values at the Gauss points.
##@item @var{options} additional options, given as pairs name/value.
##Currently only real or complex coefficient problems can be selected by
##@var{"TYPE"} and the possible values are
##@itemize
##@item @var{"real"}: all coefficients are real (default)
##@item @var{"complex"}: some coefficients might be complex
##@end itemize
##@end itemize
##
##return value
##@itemize
##@item @var{u} is the vector with the values of the solution at the nodes
##@end itemize
##
## @c Will be cut out in ??? info file and replaced with the same
## @c references explicitly there, since references to core Octave
## @c functions are not automatically transformed from here to there.
## @c BEGIN_CUT_TEXINFO
## @seealso{BVP2Dsym, BVP2eig, IBVP2D, CreateMeshRect, CreateMeshTriangle}
## @c END_CUT_TEXINFO
## @end deftypefn

function u = BVP2D(Mesh,a,b0,bx,by,f,gD,gN1,gN2,varargin)
if nargin < 9
  print_usage();
endif

Type = 'REAL'; %% default value
if (~isempty(varargin))
  for cc = 1:2:length(varargin)
    switch toupper(varargin{cc})
      case {'TYPE'}
	Type = toupper(varargin{cc+1});
      otherwise
	error('Invalid optional argument, %s. Possible value TYPE',varargin{cc});
    endswitch % switch
  endfor % for
endif % if isempty

if ((strcmp(Type,'REAL')==0)&&(strcmp(Type,'COMPLEX')==0))
  error('wrong TYPE, possible values are REAL or COMPLEX')
endif

switch Mesh.type
  case 'linear' %% first order elements
    switch Type
     case 'REAL'
       if exist("FEMEquation.oct")==3
	 [A,b] = FEMEquation(Mesh,a,b0,bx,by,f,gD,gN1,gN2);
       else
	 [A,b] = FEMEquationM(Mesh,a,b0,bx,by,f,gD,gN1,gN2);
       endif
     case 'COMPLEX'
       if exist("FEMEquationComplex.oct")==3
	 [A,b] = FEMEquationComplex(Mesh,a,b0,bx,by,f,gD,gN1,gN2);
       else
	 [A,b] = FEMEquationM(Mesh,a,b0,bx,by,f,gD,gN1,gN2);
       endif
    endswitch
  case 'quadratic'  %% second order elements
    switch Type
     case 'REAL'
       if exist("FEMEquationQuad.oct")==3
	 [A,b] = FEMEquationQuad(Mesh,a,b0,bx,by,f,gD,gN1,gN2);
       else
	 [A,b] = FEMEquationQuadM(Mesh,a,b0,bx,by,f,gD,gN1,gN2);
       endif
     case 'COMPLEX'
       if exist("FEMEquationQuadComplex.oct")==3
	 [A,b] = FEMEquationQuadComplex(Mesh,a,b0,bx,by,f,gD,gN1,gN2);
       else
	 [A,b] = FEMEquationQuadM(Mesh,a,b0,bx,by,f,gD,gN1,gN2);
       endif
    endswitch
  case 'cubic'  %% third order elements
    switch Type
     case 'REAL'
       if exist("FEMEquationCubic.oct")==3
	 [A,b] = FEMEquationCubic(Mesh,a,b0,bx,by,f,gD,gN1,gN2);
       else
	 [A,b] = FEMEquationCubicM(Mesh,a,b0,bx,by,f,gD,gN1,gN2);
       endif
     case 'COMPLEX'
       if exist("FEMEquationCubicComplex.oct")==3
	 [A,b] = FEMEquationCubicComplex(Mesh,a,b0,bx,by,f,gD,gN1,gN2);
       else
	 [A,b] = FEMEquationCubicM(Mesh,a,b0,bx,by,f,gD,gN1,gN2);
       endif
    endswitch
endswitch

u = FEMSolve(Mesh,A,b,gD);  %% solve the resulting system
endfunction
