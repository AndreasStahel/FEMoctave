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
## @deftypefn{function file}{}[@var{xnew},@var{W1mat},@var{W2mat}] = GenerateWeight1D(@var{x},@var{w1},@var{w2})
##
##   generate the weight matrices @var{W1mat} and @var{W2mat}
####
##parameters:
##@itemize
##@item @var{w1} constant, vector or function handle to evaluate the coefficient w1(x), vectorized
##@item @var{w2} (optional) constant, vector or function handle to evaluate the coefficient w2(x), vectorized
##@end itemize
##
##return values
##@itemize
##@item @var{xnew} vector with the grid points
##@item @var{W1mat} weight matrix discretizing @var{w1}
##@item @var{W2mat} (optional) weight matrix discretizing @var{w2}
##@end itemize
## @end deftypefn

function [xnew,W1mat,W2mat] = GenerateWeight1D(x,w1,w2)
if nargin<=2
  w2 = 1;
endif
  
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
w1_values = convert2values(w1);
w2_values = convert2values(w2);


W1mat = sparse(2*n+1,2*n+1); W2mat = W1mat;
%% interpolation matrix for the function values
s06 = sqrt(0.6);
G0 = [0.3+s06/2, 0.4, 0.3-s06/2;
              0,   1,         0;
      0.3-s06/2, 0.4, 0.3+s06/2];

W = [5;8;5]/18;  %% Gauss integration weights
for ind = 1:n    %% loop over all elements
  ra = 2*ind-1:2*ind+1;      %% indices for the global matrices
  ra_values = 3*ind-2:3*ind; %% indices for the values at Gauss points
  W1_elem = dx(ind)*G0'*diag(W.*w1_values(ra_values))*G0;
  W2_elem = dx(ind)*G0'*diag(W.*w2_values(ra_values))*G0;
  W1mat(ra,ra) = W1mat(ra,ra) + W1_elem;
  W2mat(ra,ra) = W2mat(ra,ra) + W2_elem;
end%for ind
end%function
