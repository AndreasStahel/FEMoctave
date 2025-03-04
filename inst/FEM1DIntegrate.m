function res = FEM1DIntegrate(x,u)
  %% Result = FEM1DIntegrate(x,u)
  %% numerical integration of FEM expression, using Simpson's rule
  %% the grid x requires the shape of FEMoctave 1D grids
  dd = diff(x(:));  w = ([dd;0]+[0;dd])/3;  w(2:2:end) *=2;
  res = w'*u(:);
endfunction
