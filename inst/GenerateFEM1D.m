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
## @deftypefn{function file}{}[@var{A},@var{M},@var{xnew}] = GenerateFEM1D(@var{x},@var{a},@var{b},@var{c},@var{d})
##
##   generate the matrices @var{A} and @var{M} to discretize the expression
##
##     -d/dx (a(x) d/dx u(x)) + b(x)*d/dx u(x) + c(x)*u(x) + d(x)*f(x)
##
##    by    A*u - M*f
##
##parameters:
##@itemize
##@item @var{a} function handle to evaluate the coefficient a(x), vectorized
##@item @var{b} function handle to evaluate the coefficient b(x), vectorized
##@item @var{c} function handle to evaluate the coefficient c(x), vectorized
##@item @var{d} function handle to evaluate the coefficient d(x), vectorized
##@end itemize
##
##return values
##@itemize
##@item @var{A} matrix discretizing the expressions involving u(x) 
##@item @var{M} matrix discretizing the evalaution of d(x)*f(x)
##@item @var{xnew} vector with the grid points
##@end itemize
##
## for an elementary demo use "demo GenerateFEM1D"
## @end deftypefn


function [A,M,xnew] = GenerateFEM1D(x,a,b,c,d)
% type  "help GenerateFEM1D"  or "demo GenerateFEM1D"

dx = diff(x(:));
n  = length(x)-1;  %% number of elements
xnew = [x(:)';x(:)' + [dx',0]/2]; xnew = xnew(:); xnew = xnew(1:end-1);
%% evaluate the coefficients at the Gauss points
xGauss = FEM1DGaussPoints(xnew);

function fun_values = convert2values(fun)  %% evaluate at Gauss points
  if (length(fun)>1)&&isnumeric(fun)  %% a given as vector
    fun_values = fun;
  elseif isnumeric(fun)               %% a given as scalar
    fun_values = fun*(ones(size(xGauss)));
  else                                %% a given as function handle
    fun_values = fun(xGauss);
  endif
endfunction
a_values = convert2values(a);
b_values = convert2values(b);
c_values = convert2values(c);
d_values = convert2values(d);


A = sparse(2*n+1,2*n+1); M = A;
%% interpolation matrix for the function values
s06 = sqrt(0.6);
G0 = [0.3+s06/2, 0.4, 0.3-s06/2;
              0,   1,         0;
      0.3-s06/2, 0.4, 0.3+s06/2];
%% interpolation matrix for the derivative values
G1 = [-1-2*s06, +4*s06, +1-2*s06;
            -1,      0,        1;
      -1+2*s06, -4*s06, +1+2*s06];

W = [5;8;5]/18;  %% Gauss integration weights
for ind = 1:n    %% loop over all elements
  ra = 2*ind-1:2*ind+1;      %% indices for the global matrices
  ra_values = 3*ind-2:3*ind; %% indices for the values at Gauss points
  M_elem = dx(ind)*G0'*diag(W.*d_values(ra_values))*G0;
  A_elem = +(G1'* diag(W.*a_values(ra_values))*G1)/dx(ind)...
           + G0'*diag(W.*b_values(ra_values))*G1 ...
	   + dx(ind)*G0'*diag(W.*c_values(ra_values))*G0;
  M(ra,ra) = M(ra,ra) + M_elem;
  A(ra,ra) = A(ra,ra) + A_elem;
end%for ind
end%function

%!demo
%! x = linspace(0,1,11);
%! a = 1; b = 0; c = 0; d = 1;  %% solve -u"=x^2  
%! [A,M,r] = GenerateFEM1D(x,a,b,c,d);
%! A = A(2:end-1,2:end-1); M = M(2:end-1,:); %% Dirichlet BC on both ends
%! f = r.^2;
%! u = A\(M*f);
%! plot(r,[0;u;0])
