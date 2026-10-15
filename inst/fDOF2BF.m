function DOF2BF = fDOF2BF(Tri)
  %% DOF2BF = fDOF2BF(Tri)
  %%
  %% determine the matrix DOF2BF to convert the DOF (Degrees Of Freedom)
  %% of a Hermite elements to the coefficients for the BF (Basis Functions)
  x = Tri(:,1);   y = Tri(:,2);  mx = sum(x)/3; my = sum(y)/3;
  x = x-mx; y = y - my;
 BF2DOF = ...
    [1 x(1) y(1) x(1)^2 x(1)*y(1) y(1)^2 x(1)^3 x(1)^2*y(1) x(1)*y(1)^2 y(1)^3;
     1 x(2) y(2) x(2)^2 x(2)*y(2) y(2)^2 x(2)^3 x(2)^2*y(2) x(2)*y(2)^2 y(2)^3;
     1 x(3) y(3) x(3)^2 x(3)*y(3) y(3)^2 x(3)^3 x(3)^2*y(3) x(3)*y(3)^2 y(3)^3;
     0  1    0   2*x(1)   y(1)      0   3*x(1)^2 2*x(1)*y(1)    y(1)^2     0;
     0  0    1     0      x(1)    2*y(1)   0      x(1)^2   2*x(1)*y(1) 3*y(1)^2;
     0  1    0   2*x(2)   y(2)      0   3*x(2)^2 2*x(2)*y(2)    y(2)^2     0;
     0  0    1     0      x(2)    2*y(2)   0      x(2)^2   2*x(2)*y(2) 3*y(2)^2;
     0  1    0   2*x(3)   y(3)      0   3*x(3)^2 2*x(3)*y(3)    y(3)^2     0;
     0  0    1     0      x(3)    2*y(3)   0      x(3)^2   2*x(3)*y(3) 3*y(3)^2;
     1  0    0     0       0       0       0        0          0          0];
  DOF2BF = inv(BF2DOF);
endfunction
