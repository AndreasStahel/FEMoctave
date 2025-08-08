function [eigval,evec,errorbound] = eigSmall(A,B,V,tol,Mode)
%  [Lambda,{Ev,err}] = eigSmall(A,V,tol,Mode)
%        solve A*Ev = Ev*diag(Lambda) standard eigenvalue problem
%
%  [Lambda,{Ev,err}] = eigSmall(A,B,V,tol,Mode)
%        solve A*Ev = B*Ev*diag(Lambda) generalized eigenvalue problem
%
%   A    is a (sparse) mxm matrix
%   B    is a (sparse) mxm matrix
%   V    is a mxn matrix, where n is the number of eigenvalues desired
%        it contains the initial eigenvectors for the iteration
%   tol  is the relative error, used as the stopping criterion
%   Mode is use to select the eigenvalues
%
%   X    is a column vector with the eigenvalues
%   EV   is a matrix whose columns represent normalized eigenvectors
%   err  is a vector with the aposteriori error estimates for the eigenvalues
%
%   this implementation is based on using eigs()
clear opts

switch nargin()
  case 1
    B = speye(size(A,1));
    opts.tol = 1e-5;
    V = 5;
  case 2
    V = B;  %% second argument is the number of eigenvalues
    B = speye(size(A,1));  % default id
    opts.tol = 1e-5;
  case 3
    opts.tol = 1e-5;
  case {4,5}
    opts.tol = tol;
endswitch

switch nargout()
  case {0,1}
    eigval = eigs(A,B,V,Mode,opts);
  case 2
    [evec,eigval] = eigs(A,B,V,Mode,opts);
    eigval = diag(eigval);
  case 3
  [evec,eigval] = eigs(A,B,V,Mode,opts);
  eigval = diag(eigval);
  errorbound = zeros(length(eigval),2);
  R = A*evec-B*evec*diag(eigval);
  R2 = B\R;
  if ~isscalar(V) V = size(V,2); endif
  for i = 1:V; %%length(errorbound)
    cc = R(:,i)'*R2(:,i);
    errorbound(i,1) = sqrt(cc);
    gap = min(abs(eigval(i)-eigval([1:i-1,i+1:end])));
    errorbound(i,2) = cc/gap;
  endfor
endswitch

%% assure that the eigenvalues are sorted
[eigval,perm] = sort(eigval);
if (nargout>=2)
  evec = evec(:,perm);
  if (nargout>=3)
    errorbound = errorbound(perm,:);
  endif
endif
