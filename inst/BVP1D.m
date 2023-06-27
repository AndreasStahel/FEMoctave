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
## @deftypefn{function file}{}[@var{x},@var{u}] = BVP1D(@var{interval},@var{a},@var{b},@var{c},@var{d},@var{f},@var{BCleft},@var{BCright})
##
##   solve a 1D boundary value problem (BVP)
##
##      -(a(x)*u'(x))' + b(x)*u'(x) + c(x)*u(x) + d(x)*f(x) = 0
##
##      with boundary conditions at the two endpoints
##@itemize
##@item Dirichlet: u(x) = g_D
##@item Neumann: a(x)*u'(x) = g_N1 + g_N2*u(x)
##@end itemize
##
##parameters:
##@itemize
##@item @var{interval} the discretized interval for the BVP
##@item @var{a} constant, vector or function handle to evaluate a(x)
##@item @var{b} constant, vector or function handle to evaluate b(x)
##@item @var{c} constant, vector or function handle to evaluate c(x)
##@item @var{d} constant, vector or function handle to evaluate d(x)
##@item @var{f} constant, vector or function handle to evaluate f(x)
##@item @var{BCleft} and @var{BCright} the two boundary conditions
##@itemize
##@item for a Dirichlet condition specify a single value @var{g_D}
##@item for a Neumann condition specify the values @var{[g_N1,g_N2]}
##@end itemize
##@end itemize
##
##return values
##@itemize
##@item @var{x} the nodes in the given interval
##@item @var{u} the values of the solution at the nodes
##@end itemize
##
## @end deftypefn

function  [x,u] = BVP1D(interval,a,b,c,d,f,BCleft,BCright)

  [A,M,x] = GenerateFEM1D(interval,a,b,c,d);

  if (length(f)>1)&&isnumeric(f)  %% f given as vector
    f_values = f;
  elseif isnumeric(f)             %% f given as scalar
    f_values = f*(ones(size(x)));
  else                            %% f given as function handle
    f_values = f(x);
  endif

  if length(BCleft)*length(BCright)==1  %% DD: Dirichlet at both ends
    a_left  = A(2:end-1,1);   %% first column
    a_right = A(2:end-1,end); %% last column
    A = A(2:end-1,2:end-1); M = M(2:end-1,:);
    u = A\(M*f_values - BCleft*a_left - BCright*a_right);
    u = [BCleft;u;BCright];
  elseif (length(BCleft)>1)&&(length(BCright)==1) %% ND: Neumann on the left, Dirichlet on the right
    a_right = A(1:end-1,end); %% last column
    A = A(1:end-1,1:end-1); M = M(1:end-1,:);
    A(1,1) += BCleft(2);
    RHS = M*f_values - BCright*a_right; RHS(1) -= BCleft(1);
    u = A\RHS;
    u = [u;BCright];
  elseif (length(BCleft)==1)&&(length(BCright)>1) %% DN: Dirichlet on the left, Neumann on the right
    a_left = A(2:end,1); %% first column
    A = A(2:end,2:end); M = M(2:end,:);
    A(end,end) -= BCright(2);
    RHS = M*f_values - BCleft*a_left; RHS(end) += BCright(1);
    u = A\RHS;
    u = [BCleft;u];
  else  %% NN: Neumann on both endpoints
    A(1,1) += BCleft(2); A(end,end) -= BCright(2);
    RHS = M*f_values; RHS(1) -= BCleft(1); RHS(end) += BCright(1);
    u = A\RHS;
  endif
