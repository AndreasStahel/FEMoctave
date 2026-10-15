function [GP,uGPmat,duxGPmat,duyGPmat] = fBF2uGPmat(Tri);
  %% [GP,uGPmat,duxGPmat,duyGPmat] = fBF2uGPmat(Tri)
  %%
  %% determine the coordinates of the Gauss points of a Hermite element
  %% and the matrices to determine u and its partial derivatives
  x = Tri(:,1); y = Tri(:,2);
  mx = mean(x); my = mean(y);
  v1 = [(x(2)+x(3))/2-x(1) (y(2)+y(3))/2-y(1)];
  v2 = [(x(3)+x(1))/2-x(2) (y(3)+y(1))/2-y(2)];
  v3 = [(x(1)+x(2))/2-x(3) (y(1)+y(2))/2-y(3)];

  l1 = (12 - 2*sqrt(15))/21;   l2 = (12 + 2*sqrt(15))/21;
  %% determine the Gauss points
  GP = [[x(1) y(1)] + l1*v1; [x(2) y(2)] + l1*v2; [x(3) y(3)] + l1*v3;
	[x(1) y(1)] + l2*v1; [x(2) y(2)] + l2*v2; [x(3) y(3)] + l2*v3;
         sum(x)/3 sum(y)/3 ];

  if nargout>1  %% evaluate the basis functions at the Gauss points
    x = GP(:,1)-mx; y = GP(:,2)-my;
    uGPmat = [ones(7,1),x,y,x.^2,x.*y,y.^2,x.^3,x.^2.*y,x.*y.^2,y.^3];
  endif
  if nargout>2  %% evaluate the first derivatives at the Gauss points
    duxGPmat = zeros(7,10); duyGPmat = zeros(7,10);
    duxGPmat(:,[2 4 5 7 8 9])  = [ones(7,1),2*x,y,3*x.^2,2*x.*y,y.^2];
    duyGPmat(:,[3 5 6 8 9 10]) = [ones(7,1),x,2*y,x.^2,2*x.*y,3*y.^2];
  endif
endfunction
