function [yi, yi_x, yi_xx] = pwquadinterp(xdata,ydata,xi)

% yi = pwquadinterp(xdata,ydata,xi)                 % to evaluate the values of the function
% [yi, yi_x, yi_xx] = pwquadinterp(xdata,ydata,xi)  % to evaluate first and second derivatives too
%
% Use the data (xdata,ydata) to determine a piecewise quadratic function
% and then evaluate this function at the points xi.
% With multiple return arguments derivatives are evaluated.
% The function requires that:
%    there are an odd number of data points,
%    xdata and xi are both in increasing order,
%    the xi values lie between xdata(1) and xdata(end).

%% original author unknown
%% modified by Andreas Stahel on May 30, 2023 to evaluate derivatives

yi = nan(size(xi));  % initialize
yi_x = yi;  yi_xx = yi;
N1 = length(xdata);
if mod(N1,2)==0
  error('Must have an odd number of data points')
end

if ~(issorted(xi) & issorted(xdata))
  error('This function assumes xi and xdata are increasing')
end

if ((xi(1)<xdata(1)) || (xi(end)>xdata(end)))
  error('This function assumes xi values lie within the interval of xdata')
end

ii = 1;
for j=3:2:N1
   % determine polynomial on interval [xdata(j-2), xdata(j)]
   x1 = xdata(j-2);  x2 = xdata(j-1);   x3 = xdata(j);
   f1 = ydata(j-2);  f2 = ydata(j-1);   f3 = ydata(j);
   % compute divided differences:
   f12 = (f2-f1) / (x2-x1);
   f23 = (f3-f2) / (x3-x2);
   f123 = (f23-f12) / (x3-x1);

% set function values and derivatives for any xi(ii) that lie in this interval:
   while ((ii<=length(xi)) && (xi(ii)<=x3) )
      yi(ii) = f1 + f12*(xi(ii)-x1) + f123*(xi(ii)-x1)*(xi(ii)-x2);
      yi_x(ii) = f12 + f123*((xi(ii)-x1)+(xi(ii)-x2));
      yi_xx(ii) = 2*f123;
      ii = ii+1;
   end%while
end%for

if ii<=length(xi)   % we failed to set all the yi values
   disp('*** This should not happen due to checks on input!')
end%if


