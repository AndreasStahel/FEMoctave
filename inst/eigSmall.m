function [eigval,evec,errorbound] = eigSmall(A,B,V,tol)
%  [Lambda,{Ev,err}] = eigSmall(A,V,tol)     
%        solve A*Ev = Ev*diag(Lambda) standard eigenvalue problem
%
%  [Lambda,{Ev,err}] = eigSmall(A,B,V,tol)   
%        solve A*Ev = B*Ev*diag(Lambda) generalized eigenvalue problem
%
%   A   is mxt, where t-1 is number of non-zero superdiagonals
%   B   is mxs, where s-1 is number of non-zero superdiagonals
%   V   is mxn, where n is the number of eigenvalues desired
%       contains the initial eigenvectors for the iteration
%   tol is the relative error, used as the stopping criterion
%
%   X   is a column vector with the eigenvalues
%   EV  is a matrix whose columns represent normalized eigenvectors
%   err is a vector with the aposteriori error estimates for the eigenvalues
%
%   this implementation is based on using eigs()

clear opts
if nargin == 4
  opts.tol=tol;
else
  opts.dummy=1;
endif

if nargout() == 1
  eigval = eigs(A,B,V,'sm',opts);
endif

if nargout() == 2
  [evec,eigval] = eigs(A,B,V,'sm',opts);
  eigval = diag(eigval)
endif

if nargout() == 3
  [evec,eigval] = eigs(A,B,V,'sm',opts);
  eigval = diag(eigval);
  errorbound = zeros(length(eigval),2);
  R = A*evec-B*evec*diag(eigval);
  R2 = B\R;
  for i = 1:V; %%length(errorbound)
    cc = R(:,i)'*R2(:,i);
    errorbound(i,1) = sqrt(cc);
    gap = min(abs(eigval(i)-eigval([1:i-1,i+1:end])));
    errorbound(i,2) = cc/gap;
  endfor
endif

%% assure that the eigenvalues are sorted
[eigval,perm] = sort(eigval);
if (nargout>=2)
  evec = evec(:,perm);
  if (nargout>=3)
    errorbound = errorbound(perm,:);
  endif
endif


