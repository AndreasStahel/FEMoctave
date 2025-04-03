function wMat = GenerateWeight2D(Mesh)
%[...] = GenerateWeight2D (...)
%
%  evalutes the matrix W required in BVP2DNL() to solve
%  nonlinear boundary value problems in 2D 
%
%  W = GenerateWeight2D(Mesh)
%
%    Mesh is the mesh describing the domain
%
%    W   is the matrix to be used by BVP2DNL()
%

nElem = size(Mesh.elem,1); nGP = size(Mesh.GP,1)  ;

%% create memory for the sparse matrix and the RHS vector

switch Mesh.type
  case 'linear'
    Wi = zeros(9*nElem,1); Wj = Wi; Wval = Wi;%% maximal number of contributions
    M = [4 1 1; 1 4 1; 1 1 4]/6;  %% interpolation matrix, is symmetric
      % insert the element matrices and vectors into the global matrix
    Counter = 1;  %% counter for the entries in the matrix 
    for k = 1:nElem   %%for each element
      area = Mesh.elemArea(k);  % area = 0.5*det(T)
      El_mat = area/3*M;
      dofs = Mesh.node2DOF(Mesh.elem(k,:));
      for k1 = 1:3
	if dofs(k1)>0 % k1 is free node
	  Wi(Counter:Counter+2)   = dofs(k1);
	  Wj(Counter:Counter+2)   = (k-1)*3+[1:3];
	  Wval(Counter:Counter+2) = El_mat(k1,:);
	  Counter = Counter+3;
	endif  % dofs(k1)
      endfor % k1
    endfor % k (elements)    
    %% add up to create the sparse matrix
    Wi = Wi(1:Counter-1); Wj = Wj(1:Counter-1); Wval = Wval(1:Counter-1);
    wMat = sparse(Wi,Wj,Wval,Mesh.nDOF,nGP);

  case 'quadratic'
    Wi = zeros(42*nElem,1); Wj=Wi; Wval=Wi; %% maximal number of contributions
    l1 = (12 - 2*sqrt(15))/21;   l2 = (12 + 2*sqrt(15))/21;
    w1 = (155 - sqrt(15))/2400;  w2 = (155 + sqrt(15))/2400;
    w3 = 0.1125;        w = [w1,w1,w1,w2,w2,w2,w3]';
    xi = [l1/2, 1-l1, l1/2, l2/2, 1-l2, l2/2, 1/3]';
    nu = [l1/2, l1/2, 1-l1, l2/2, l2/2, 1-l2, 1/3]';
    %% the interpolation matrix for the function
    M = [(1-xi-nu).*(1-2*xi-2*nu) xi.*(2*xi-1) nu.*(2*nu-1) 4*xi.*nu 4*nu.*(1-xi-nu) 4*xi.*(1-xi-nu)]';
    %% insert the element matrices and vectors into the global matrix
    Counter = 1;  %% counter for the entries in the matrix 
    for k = 1:nElem   %%for each element
      area = Mesh.elemArea(k);  % area = 0.5*det(T)
      El_mat = 2*area*M*diag(w);
      dofs = Mesh.node2DOF(Mesh.elem(k,:));
      for k1 = 1:6
	if dofs(k1)>0 % k1 is free node
	  Wi(Counter:Counter+6)   = dofs(k1);
	  Wj(Counter:Counter+6)   = (k-1)*7+[1:7];
	  Wval(Counter:Counter+6) = El_mat(k1,:);
	  Counter = Counter+7;
	endif  % dofs(k1)
      endfor % k1
    endfor % k (elements)    
    %% add up to create the sparse matrix
    Wi = Wi(1:Counter-1); Wj = Wj(1:Counter-1); Wval = Wval(1:Counter-1);
    wMat = sparse(Wi,Wj,Wval,Mesh.nDOF,nGP);

  case 'cubic'
    Wi = zeros(70*nElem,1); Wj=Wi; Wval=Wi; %% maximal number of contributions
    l1 = (12 - 2*sqrt(15))/21;   l2 = (12 + 2*sqrt(15))/21;
    w1 = (155 - sqrt(15))/2400;  w2 = (155 + sqrt(15))/2400;
    w3 = 0.1125;        w = [w1,w1,w1,w2,w2,w2,w3]';
    xi = [l1/2, 1-l1, l1/2, l2/2, 1-l2, l2/2, 1/3]';
    nu = [l1/2, l1/2, 1-l1, l2/2, l2/2, 1-l2, 1/3]';
    %% the interpolation matrix for the function
    M = [1-11/2*xi-11/2*nu+9*xi.^2+18*xi.*nu+9*nu.^2-9/2*xi.^3-27/2*xi.^2.*nu-27/2*xi.*nu.^2-9/2*nu.^3,...
     xi-9/2*xi.^2+9/2*xi.^3,...
     nu-9/2*nu.^2+9/2*nu.^3,...
     -9/2*xi.*nu+27/2*xi.^2.*nu,...
     -9/2*xi.*nu+27/2*xi.*nu.^2,...
     -9/2*nu+9/2*xi.*nu+18*nu.^2-27/2*xi.*nu.^2-27/2*nu.^3,...
     9*nu-45/2*xi.*nu-45/2*nu.^2+27/2*xi.^2.*nu+27*xi.*nu.^2+27/2*nu.^3,...
     9*xi-45/2*xi.^2-45/2*xi.*nu+27/2*xi.^3+27*xi.^2.*nu+27/2*xi.*nu.^2,...
     -9/2*xi+18*xi.^2+9/2*xi.*nu-27/2*xi.^3-27/2*xi.^2.*nu,...
     27*xi.*nu-27*xi.^2.*nu-27*xi.*nu.^2 ]';
    %% insert the element matrices and vectors into the global matrix
    Counter = 1;  %% counter for the entries in the matrix 
    for k = 1:nElem   %%for each element
      area = Mesh.elemArea(k);  % area = 0.5*det(T)
      El_mat = 2*area*M*diag(w);
      dofs = Mesh.node2DOF(Mesh.elem(k,:));
      for k1 = 1:10
	if dofs(k1)>0 % k1 is free node
	  Wi(Counter:Counter+6)   = dofs(k1);
	  Wj(Counter:Counter+6)   = (k-1)*7+[1:7];
	  Wval(Counter:Counter+6) = El_mat(k1,:);
	  Counter = Counter+7;
	endif  % dofs(k1)
      endfor % k1
    endfor % k (elements)    
    %% add up to create the sparse matrix
    Wi = Wi(1:Counter-1); Wj = Wj(1:Counter-1); Wval = Wval(1:Counter-1);
    wMat = sparse(Wi,Wj,Wval,Mesh.nDOF,nGP);
endswitch

    
endfunction
