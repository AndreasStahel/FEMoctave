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
## Created: 2023-05-30


## -*- texinfo -*-
## @deftypefn{function file}{}[@var{x},@var{Eval},@var{Evec},@var{errorbound}] = BVP1Deig(@var{interval},@var{a},@var{b},@var{c},@var{w},@var{BCleft},@var{BCright},@var{nVec},@var{tol})
##
##determine the smallest eigenvalues @var{Eval} and eigenfunctions @var{Evec} for the BVP
##
##      -(a(x)*u'(x))' + b(x)*u'(x) + c(x)*u(x) = eVal*w(x)*u(x)
##
##                                           u  = 0      on Dirichlet boundary
##
##                                        a*u'  = g_N2*u on  Neumann boundary
##
##parameters:
##@itemize
##@item @var{interval} the discretized interval for the BVP
##@item @var{a} constant, vector or function handle to evaluate a(x)
##@item @var{b} constant, vector or function handle to evaluate b(x)
##@item @var{c} constant, vector or function handle to evaluate c(x)
##@item @var{w} constant, vector or function handle to evaluate w(x)
##@item @var{BCleft} and @var{BCright} the two boundary conditions
##@itemize
##@item for a Dirichlet condition specify a single value @var{g_D=0}
##@item for a Neumann condition specify the values @var{[g_N1,g_N2]=[0,g_N2]}
##@end itemize
##@item @var{nVec} the number of smallest eigenvalues to be computed
##@item@var{tol} tolerance for the eigenvalue iteration, default 1e-5
##@end itemize
##
##return values
##@itemize
##@item @var{x} the nodes in the given interval
##@item @var{eVal} the eigenvalues
##@item @var{eVec} the matrix of eigenvectors of the solutions at the nodes
##@item @var{errorbound} a matrix with error bounds of the eigenvalues
##@end itemize
##
## @end deftypefn

function  [x,eVal,eVec,errorbound] = BVP1Deig(interval,a,b,c,w,BCleft,BCright,nVec,tol)

if (nargin==8) tol = 1e-5;endif

[A,M,x] = GenerateFEM1D(interval,a,b,c,w);

if length(BCleft)*length(BCright)==1  %% DD: Dirichlet at both ends
  A = A(2:end-1,2:end-1); M = M(2:end-1,2:end-1);
  if (nargout>=2)
    if (nargout==4)
      [eVal,eVec,errorbound] = eigSmall(A,M,nVec,tol);
    else
      [eVal,eVec]            = eigSmall(A,M,nVec,tol);
    endif
  endif %% nargout>3
  eVec = [zeros(1,nVec);eVec;zeros(1,nVec)];
elseif (length(BCleft)>1)&&(length(BCright)==1) %% ND: Neumann on the left, Dirichlet on the right
    A = A(1:end-1,1:end-1); M = M(1:end-1,1:end-1);
    A(1,1) += BCleft(2);
    if (nargout>=3)
      if (nargout==4)
	[eVal,eVec,errorbound] = eigSmall(A,M,nVec,tol);
      else
	[eVal,eVec]            = eigSmall(A,M,nVec,tol);
      endif
    endif %% nargout>3
    eVec = [eVec;zeros(1,nVec)];
elseif (length(BCleft)==1)&&(length(BCright)>1) %% DN: Dirichlet on the left, Neumann on the right
  A = A(2:end,2:end); M = M(2:end,2:end);
  A(end,end) -= BCright(2);
  if (nargout>=3)
    if (nargout==4)
      [eVal,eVec,errorbound] = eigSmall(A,M,nVec,tol);
    else
      [eVal,eVec]            = eigSmall(A,M,nVec,tol);
    endif
  endif %% nargout>3
  eVec = [zeros(1,nVec);eVec];
else  %% NN: Neumann on both endpoints
  A(1,1) += BCleft(2); A(end,end) -= BCright(2);
  pi
  if (nargout>=3)
    if (nargout==4)
      [eVal,eVec,errorbound] = eigSmall(A,M,nVec,tol);
    else
      [eVal,eVec]            = eigSmall(A,M,nVec,tol);
    endif
  endif %% nargout>3
endif
