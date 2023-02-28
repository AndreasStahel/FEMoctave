function [gMat,gVec] = FEMEquation(Mesh,aFunc,b0Func,bxFunc,byFunc,fFunc,gDFunc,gN1Func,gN2Func)
%[...] = FEMEquation (...)
%  sets up the system of linear equations for a numerical solution of a PDE
%
%  [A,b,n2d] = FEMEquation(Mesh,'a','b0','bx','by','f','gD','gN1','gN2')
%    Mesh is the mesh describing the domain\n\
%         see ReadMesh() for the description of the format
%   'a','b0','bx','by','f','gD','gN1','gN2' are the functions and coefficients 
%         for the boundary value problem. They can be given as a scalar value
%         or as a string with the function name 
%
%  -div(a*grad u-u*(bx,by)) + b0*u = f    in domain
%                          u = gD         on Dirichlet section of the boundary
%         (a*du-u*(bx,by))*n = gN1+gN2*u  on Neumann section of the boundary 
%
%
% A   is the matrix of the system to be solved.
% b   is the RHS of the system to be solved.
%

%% evaluate the functions a, b and f

nElem = size(Mesh.elem,1); nGP   = size(Mesh.GP,1)  ;

if ischar(aFunc)
  aV = reshape(feval(aFunc,Mesh.GP,Mesh.GPT),nGP/nElem,nElem);
elseif isscalar(aFunc)
  aV = aFunc*ones(nGP/nElem,nElem);
else
  aV = reshape(aFunc,nGP/nElem,nElem);
endif

if ischar(b0Func)
  bV = reshape(feval(b0Func,Mesh.GP,Mesh.GPT),nGP/nElem,nElem);
elseif isscalar(b0Func)
  bV = b0Func*ones(nGP/nElem,nElem);
else
  bV = reshape(b0Func,nGP/nElem,nElem);
endif

if ischar(bxFunc)
  bxV = reshape(feval(bxFunc,Mesh.GP,Mesh.GPT),nGP/nElem,nElem);
elseif isscalar(bxFunc)
  bxV = bxFunc*ones(nGP/nElem,nElem);
else
  bxV = reshape(bxFunc,nGP/nElem,nElem);
endif

if ischar(byFunc)
  byV = reshape(feval(byFunc,Mesh.GP,Mesh.GPT),nGP/nElem,nElem);
elseif isscalar(byFunc)
  byV = byFunc*ones(nGP/nElem,nElem);
else
  byV = reshape(byFunc,nGP/nElem,nElem);
endif

if ischar(fFunc)
  fV = reshape(feval(fFunc,Mesh.GP,Mesh.GPT),nGP/nElem,nElem);
elseif isscalar(fFunc)
  fV = fFunc*ones(nGP/nElem,nElem);
else
  fV = reshape(fFunc,nGP/nElem,nElem);
endif

%% create memory for the sparse matrix and the RHS vector

Si   = zeros(9*nElem,1); Sj = Si; Sval = Si;  %% maximal number of contributions
gVec = zeros(Mesh.nDOF,1);

M = [4 1 1; 1 4 1; 1 1 4]/6;  %% interpolation matrix, is symmetric
% insert the element matrices and vectors into the global matrix
ptrDOF = 1;  %% counter for the DOF we are working on
for k = 1:nElem   %%for each element
  cor = Mesh.nodes(Mesh.elem(k,:),:);  % coordinates of the nodes
 %% compute element stiffness matrix and vector
  area = Mesh.elemArea(k);  % area = 0.5*det(T)
  G = [cor(3,2)-cor(2,2),cor(1,2)-cor(3,2),cor(2,2)-cor(1,2);...
       cor(2,1)-cor(3,1),cor(3,1)-cor(1,1),cor(1,1)-cor(2,1)];
  Mdiffusion = (bxV(:,k)*G(1,:) + byV(:,k)*G(2,:))'*M/6;
  mat = sum(aV(:,k))/(12*area)*G'*G + area/3*M*diag(bV(:,k))*M+Mdiffusion;
  vec = -area/3*M*fV(:,k);
  dofs = Mesh.node2DOF(Mesh.elem(k,:));
  for k1 = 1:3
    if dofs(k1)>0 % k1 is free node
      gVec(dofs(k1)) += vec(k1);
      for k2 = 1:3
	if dofs(k2)>0  % k2 is free node
 	  %% gMat(dofs(k1),dofs(k2)) = gMat(dofs(k1),dofs(k2)) + mat(k1,k2);
	  Si(ptrDOF)   = dofs(k1);
	  Sj(ptrDOF)   = dofs(k2);
	  Sval(ptrDOF) = mat(k1,k2);
	  ptrDOF++;
	else  %% k2 is a Dirichlet node
	  if   ischar(gDFunc) gD = feval(gDFunc,cor(k2,:));
	  else                gD = gDFunc;
	  endif% ischchar
 	  gVec(dofs(k1)) += mat(k1,k2)*gD;
	endif % dofs(k2)>0
      endfor  % k2
    endif  % dofs(k1)
  endfor % k1
endfor % k (elements)

%% add up to create the sparse matrix
%%Sval(find(abs(Sval)<1e-15)) = 0;  % eliminate very small entries
Si = Si(1:ptrDOF-1); Sj = Sj(1:ptrDOF-1); Sval = Sval(1:ptrDOF-1);
gMat = sparse(Si,Sj,Sval,Mesh.nDOF,Mesh.nDOF);

%% insert the edge contributions
for k = 1:size(Mesh.edges)(1)
  if Mesh.edgesT(k)<-1  % it is a Neumann edge
    cor = Mesh.nodes(Mesh.edges(k,:),:);
%%    edgeVec = ElementContributionEdge(corners,gNFunc);
    p1 = (cor(1,:)+cor(2,:))/2 - (cor(2,:)-cor(1,:))/(2*sqrt(3)); 
    p2 = (cor(1,:)+cor(2,:))/2 + (cor(2,:)-cor(1,:))/(2*sqrt(3)); 
    if ischar(gN1Func) g = feval(gN1Func,[p1;p2]);
    else               g = gN1Func*ones(2,1);
    endif

    alpha = (1-1/sqrt(3))/2; L = norm(cor(2,:)-cor(1,:))/2;
    edgeVec = L*[(1-alpha)*g(1)+alpha*g(2);
		 alpha*g(1)+(1-alpha)*g(2)];

    if ischar(gN2Func) g = feval(gN2Func,[p1;p2]);  %% evaluate gN2
    else               g = gN2Func*ones(2,1);
    endif
    B = L*[1-alpha alpha;alpha 1-alpha]*diag(g)*[1-alpha alpha;alpha 1-alpha];

    if ischar(gDFunc) g = feval(gDFunc,cor);  %% evaluate Dirichlet values
    else              g = gDFunc*ones(2,1);
    endif
    
    dofs = Mesh.node2DOF(Mesh.edges(k,:));
    if dofs(1)>0   %% node 1 free
      if dofs(2)>0 %% both points free
	gVec(dofs) -= edgeVec;
        gMat(dofs,dofs) -= B;
      else  % node 1 free, node 2 Dirichlet
	gVec(dofs(1)) -= edgeVec(1) + B(1,2)*g(2);
	gMat(dofs(1),dofs(1)) -= B(1,1);
      endif% dofs(2)>0
    else  %% node 1 Dirichlet, node 2 free
      gVec(dofs(2)) -= edgeVec(2) + B(2,1)*g(1); 
      gMat(dofs(2),dofs(2)) -= B(2,2);
    endif % dofs(1)>0
  endif % Neumann edge
endfor  % k edges
endfunction
