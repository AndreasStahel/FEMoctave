## Copyright (C) 2020-2025 Andreas Stahel
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
## Created: 2025-01-08

## -*- texinfo -*-
## @deftypefn{function file}{}[@var{Eval},@var{Evec},@var{errorbound}] = BVP2Deig(@var{Mesh},@var{a},@var{b0},@var{w},@var{gN2},@var{nVec},@var{tol},@var{options})
##
##determine the smallest eigenvalues @var{Eval} and eigenfunctions @var{Evec} for the BVP
##
##@verbatim
##     -div(a*grad u) + b0*u = Eval*w*u  in domain
##                        u  = 0         on Dirichlet boundary
##             n*(a*grad u)  = gN2*u on  Neumann boundary
##@end verbatim
##
##parameters:
##@itemize
##@item @var{Mesh} the mesh describing the domain and the boundary types
##@item @var{a},@var{b0},@var{w},@var{gN2}
## the coefficients and functions describing the PDE.
##@*Any constant function can be given by its scalar value.
##@*The functions @var{a},@var{b0} and @var{w} may also be given as vectors
##with the values of the function at the Gauss points.
##@item The coefficient @var{a} can also be a symmetric matrix
## A=[axx,axy;axy,ayy] given by the row vector [axx,ayy,axy].
## It can be given as row vector or as string with the function name or
## as nx3 matrix with the values at the Gauss points.
##@item@var{nVec} number of smallest eigenvalues to be computed
##@item @var{options} additional options, given as pairs name/value.
##@itemize
##@item@var{tol} tolerance for the eigenvalue iteration, default 1e-5
##@item@var{"TYPE"} select real or complex coefficients
##@itemize
##@item @var{"real"}: all coefficients are real (default)
##@item @var{"complex"}: some coefficients might be complex
##@end itemize
##@end itemize
##@end itemize
##
##return values:
##@itemize
##@item @var{Eval} is the vector with the eigenvalues
##@item @var{Evec} is the matrix with the eigenvectors as columns
##@item @var{errorbound} is a matrix with error bounds of the eigenvalues
##@end itemize
##
## @c Will be cut out in ??? info file and replaced with the same
## @c references explicitly there, since references to core Octave
## @c functions are not automatically transformed from here to there.
## @c BEGIN_CUT_TEXINFO
## @seealso{BVP2D, BVP2sym, IBVP2D, IBVP2Dsym, I2BVP2D}
## @c END_CUT_TEXINFO
## @end deftypefn

function [la,resVec,errorbound] = BVP2Deig(Mesh,aFunc,bFunc,wFunc,gN2Func,nVec,varargin)
%  [la,ev] = BVP2Deig(Mesh,a,b0,w,gN2,nVec,options)
%
%  determine the smallest eigenvalues la and eigenfunctions ev for the BVP
%
%  -div(a*grad u) + b0*u = la*w*u  in domain
%                      u = 0       on Dirichlet boundary
%           n*(a*grad u) = gN2*u   on Neumann boundary
%
%  la      = FEMEig(Mesh,aFunc,bFunc,wFunc,gN2Func,nVec,tol,options)
%  [la,ev] = FEMEig(Mesh,aFunc,bFunc,wFunc,gN2Func,nVec,tol,options)
%
% nodes elem, edges contains information about the mesh
%         see ReadMesh() for the description of the format
% a, b0, w, gN2 are function files for the coefficient functions
%         may also be given as vectors or scalar values
% nVec    is the number of eigenvalues to be computed
% tol    
%         if not given tol = 1e-5 is used as default
% options additional options, given as pairs name/value.
%    * TOL  tolerance for error of the eigenvalues, default 1e-5
%    * TYPE select real or complex coefficients
%      * "REAL": all coefficients are real (default)
%      * "COMPLEX": some coefficients might be complex
%
% la  is the vector containing the eigenvalues
% ev  is the matrix with the eigenvectors as columns
%
%see also ????

if ((nargin<6))
  help("BVP2Deig");
  print_usage()
endif

Type = 'REAL'; %% default value
tol  = 1e-5  ; %% default value
if (~isempty(varargin))
  for cc = 1:2:length(varargin)
    switch toupper(varargin{cc})
      case {'TYPE'}
	Type = toupper(varargin{cc+1});
      case {'TOL'}
	tol = varargin{cc+1};
      otherwise
	error('Invalid optional argument, %s. Possible values: TOL, TYPE',varargin{cc});
    endswitch % switch
  endfor % for
endif % if isempty

if ((strcmp(Type,'REAL')==0)&&(strcmp(Type,'COMPLEX')==0))
  error('wrong TYPE, possible values are REAL or COMPLEX')
endif

switch Mesh.type
  case 'linear'  %% first order elements
    switch Type
     case 'REAL'
       if exist("FEMEquation.oct")==3
	 [aMat,~] = FEMEquation(Mesh,aFunc,bFunc,0,0,0,0,0,gN2Func);
	 [wMat,~] = FEMEquation(Mesh,0    ,wFunc,0,0,0,0,0,0);
       else
	 [aMat,~] = FEMEquationM(Mesh,aFunc,bFunc,0,0,0,0,0,gN2Func);
	 [wMat,~] = FEMEquationM(Mesh,0    ,wFunc,0,0,0,0,0,0);
       endif
     case 'COMPLEX'
       if exist("FEMEquationComplex.oct")==3
	 [aMat,~] = FEMEquationComplex(Mesh,aFunc,bFunc,0,0,0,0,0,gN2Func);
	 [wMat,~] = FEMEquationComplex(Mesh,0    ,wFunc,0,0,0,0,0,0);
       else
	 [aMat,~] = FEMEquationM(Mesh,aFunc,bFunc,0,0,0,0,0,gN2Func);
	 [wMat,~] = FEMEquationM(Mesh,0    ,wFunc,0,0,0,0,0,0);
       endif
    endswitch
  case 'quadratic' %% second order elements
    switch Type
     case 'REAL'
       if exist("FEMEquationQuad.oct")==3
	 [aMat,~] = FEMEquationQuad(Mesh,aFunc,bFunc,0,0,0,0,0,gN2Func);
	 [wMat,~] = FEMEquationQuad(Mesh,0    ,wFunc,0,0,0,0,0,0);
       else
	 [aMat,~] = FEMEquationQuadM(Mesh,aFunc,bFunc,0,0,0,0,0,gN2Func);
	 [wMat,~] = FEMEquationQuadM(Mesh,0    ,wFunc,0,0,0,0,0,0);
       endif
     case 'COMPLEX'
       if exist("FEMEquationQuadComplex.oct")==3
	 [aMat,~] = FEMEquationQuadComplex(Mesh,aFunc,bFunc,0,0,0,0,0,gN2Func);
	 [wMat,~] = FEMEquationQuadComplex(Mesh,0    ,wFunc,0,0,0,0,0,0);
       else
	 [aMat,~] = FEMEquationQuadM(Mesh,aFunc,bFunc,0,0,0,0,0,gN2Func);
	 [wMat,~] = FEMEquationQuadM(Mesh,0    ,wFunc,0,0,0,0,0,0);
       endif
    endswitch
  case 'cubic' %% third order elements
    switch Type
     case 'REAL'
       if exist("FEMEquationCubic.oct")==3
	 [aMat,~] = FEMEquationCubic(Mesh,aFunc,bFunc,0,0,0,0,0,gN2Func);
	 [wMat,~] = FEMEquationCubic(Mesh,0    ,wFunc,0,0,0,0,0,0);
       else
	 [aMat,~] = FEMEquationCubicM(Mesh,aFunc,bFunc,0,0,0,0,0,gN2Func);
	 [wMat,~] = FEMEquationCubicM(Mesh,0    ,wFunc,0,0,0,0,0,0);
       endif
     case 'COMPLEX'
       if exist("FEMEquationCubicComplex.oct")==3
	 [aMat,~] = FEMEquationCubicComplex(Mesh,aFunc,bFunc,0,0,0,0,0,gN2Func);
	 [wMat,~] = FEMEquationCubicComplex(Mesh,0    ,wFunc,0,0,0,0,0,0);
       else
	 [aMat,~] = FEMEquationCubicM(Mesh,aFunc,bFunc,0,0,0,0,0,gN2Func);
	 [wMat,~] = FEMEquationCubicM(Mesh,0    ,wFunc,0,0,0,0,0,0);
       endif
    endswitch
endswitch

if (nargout==1)
  la = eigSmall(aMat,wMat,nVec,tol);
endif

if (nargout>=2)
  if (nargout==3)
    [la,ug,errorbound] = eigSmall(aMat,wMat,nVec,tol);
  else
    [la,ug]            = eigSmall(aMat,wMat,nVec,tol);
  endif
  n = length(Mesh.node2DOF);
  m = size(ug)(2);

  resVec = zeros(n,m);
  [Ind,~,Val] = find(Mesh.node2DOF>0);
  resVec(Ind,:) = ug(Val ,:);
endif
endfunction
