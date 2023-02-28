function [gMat,gVec] = FEMEquationCubic(Mesh,aFunc,b0Func,bxFunc,byFunc,fFunc,gDFunc,gN1Func,gN2Func)
%[...] = FEMEquationCubic (...)
%  sets up the system of linear equations for a numerical solution of a PDE
%
%  [A,b] = FEMEquationCubicM(Mesh,'a','b0','bx','by','f','gD','gN1','gN2')
%    Mesh is the mesh describing the domain\n\
%         see ReadMesh() for the description of the format
%   'a','b0','bx',by','f','gD','gN1','gN2' are the functions and coefficients 
%         for the boundary value problem. They can be given as a scalar value
%         or as a string with the function name 
%
%  -div(a*grad u-u*(bx,by)) + b0*u = f    in domain%
%                          u = gD         on Dirichlet section of the boundary
%         (a*du-u*(bx,by))*n = gN1+gN2*u  on Neumann section of the boundary 
%
% A   is the matrix of the system to be solved.
% b   is the RHS of the system to be solved.
%

%% evaluate the functions a b and f

nElem = size(Mesh.elem,1);
nGP   = size(Mesh.GP,1)  ;

if ischar(aFunc)
  aV = reshape(feval(aFunc,Mesh.GP,Mesh.nodesT),nGP/nElem,nElem);
elseif isscalar(aFunc)
  aV = aFunc*ones(nGP/nElem,nElem);
else
  aV = reshape(aFunc,nGP/nElem,nElem);
endif

if ischar(b0Func)
  bV = reshape(feval(b0Func,Mesh.GP,Mesh.nodesT),nGP/nElem,nElem);
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
  fV = reshape(feval(fFunc,Mesh.GP,Mesh.nodesT),nGP/nElem,nElem);
elseif isscalar(fFunc)
  fV = fFunc*ones(nGP/nElem,nElem);
else
  fV = reshape(fFunc,nGP/nElem,nElem);
endif

%% create memory for the sparse matrix and the RHS vector
Si  = zeros(100*nElem,1); Sj = Si; Sval = Si; %% maximal number of contributions
gVec= zeros(Mesh.nDOF,1);

l1 = (12 - 2*sqrt(15))/21;   l2 = (12 + 2*sqrt(15))/21;
w1 = (155 - sqrt(15))/2400;  w2 = (155 + sqrt(15))/2400;
w3 = 0.1125;         w = [w1,w1,w1,w2,w2,w2,w3]';

xi = [l1/2, 1-l1, l1/2, l2/2, 1-l2, l2/2, 1/3]';
nu = [l1/2, l1/2, 1-l1, l2/2, l2/2, 1-l2, 1/3]';

%% the interpolation matrices for the function and its partial derivatives
M = [1-11/2*xi-11/2*nu+9*xi.^2+18*xi.*nu+9*nu.^2-9/2*xi.^3-27/2*xi.^2.*nu-27/2*xi.*nu.^2-9/2*nu.^3,...
     xi-9/2*xi.^2+9/2*xi.^3,...
     nu-9/2*nu.^2+9/2*nu.^3,...
     -9/2*xi.*nu+27/2*xi.^2.*nu,...
     -9/2*xi.*nu+27/2*xi.*nu.^2,...
     -9/2*nu+9/2*xi.*nu+18*nu.^2-27/2*xi.*nu.^2-27/2*nu.^3,...
     9*nu-45/2*xi.*nu-45/2*nu.^2+27/2*xi.^2.*nu+27*xi.*nu.^2+27/2*nu.^3,...
     9*xi-45/2*xi.^2-45/2*xi.*nu+27/2*xi.^3+27*xi.^2.*nu+27/2*xi.*nu.^2,...
     -9/2*xi+18*xi.^2+9/2*xi.*nu-27/2*xi.^3-27/2*xi.^2.*nu,...
     27*xi.*nu-27*xi.^2.*nu-27*xi.*nu.^2 ];
z = zeros(size(xi));

Mxi = [-11/2+18*xi+18*nu-27/2*xi.^2-27*xi.*nu-27/2*nu.^2,...
       1-9*xi+27/2*xi.^2,...
       z,...
       -9/2*nu+27*xi.*nu,...
       -9/2*nu+27/2*nu.^2,...
       9/2*nu-27/2*nu.^2,...
       -45/2*nu+27*xi.*nu+27*nu.^2,...
       9-45*xi-45/2*nu+81/2*xi.^2+54*xi.*nu+27/2*nu.^2,...
       -9/2+36*xi+9/2*nu-81/2*xi.^2-27*xi.*nu,...
       27*nu-54*xi.*nu-27*nu.^2 ];

Mnu = [-11/2+18*xi+18*nu-27/2*xi.^2-27*xi.*nu-27/2*nu.^2,...
       z,...
       1-9*nu+27/2*nu.^2,...
       -9/2*(1-3*xi).*xi,...
       -9/2*xi+27*nu.*xi,...
       -9/2+9/2*xi+36*nu-27*xi.*nu-81/2*nu.^2,...
       9-45/2*xi-45*nu+27/2*xi.^2+54*xi.*nu+81/2*nu.^2,...
       -45/2*xi+27*xi.^2+27*xi.*nu,...
       +9/2*xi-27/2*xi.^2,...
       27*xi-27*xi.^2-54*xi.*nu ];

% insert the element matrices and vectors into the global matrix
ptrDOF = 1;       %% counter for the DOF we are working on
for k = 1:nElem   %%for each element
  cor = Mesh.nodes(Mesh.elem(k,:),:);  % coordinates of the nodes
 %% compute element stiffness matrix and vector
  area = Mesh.elemArea(k);  % area = 0.5*det(T)
  G = [cor(3,2)-cor(2,2),cor(1,2)-cor(3,2),cor(2,2)-cor(1,2);...
       cor(2,1)-cor(3,1),cor(3,1)-cor(1,1),cor(1,1)-cor(2,1)];
  B = M'*diag(w.*bV(:,k))*M;
  Mtemp = -G(1,2)*Mxi-G(1,3)*Mnu;
  Ax = Mtemp'*diag(w.*aV(:,k))*Mtemp;
  Mtemp = -G(2,2)*Mxi-G(2,3)*Mnu;
  Ay = Mtemp'*diag(w.*aV(:,k))*Mtemp;
  Abxy = ((G(1,2)*Mxi + G(1,3)*Mnu)'*diag(w.*bxV(:,k))...
         +(G(2,2)*Mxi + G(2,3)*Mnu)'*diag(w.*byV(:,k)) )*M;
  mat = 0.5/area*(Ax+Ay) + 2*area*B + Abxy;
  vec = -2*area*M'*(w.*fV(:,k));
  dofs = Mesh.node2DOF(Mesh.elem(k,:));
  for k1 = 1:10
    if dofs(k1)>0 % k1 is free node
      gVec(dofs(k1)) +=  vec(k1);
      for k2 = 1:10
	if (dofs(k2)>0)  % k2 is free node
 	  %% gMat(dofs(k1),dofs(k2)) = gMat(dofs(k1),dofs(k2)) + mat(k1,k2);
	  Si(ptrDOF) = dofs(k1); Sj(ptrDOF) = dofs(k2);
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
Sval(find(abs(Sval)<1e-15))=0;  % eliminate very small entries
Si = Si(1:ptrDOF-1); Sj = Sj(1:ptrDOF-1); Sval = Sval(1:ptrDOF-1);
gMat = sparse(Si,Sj,Sval,Mesh.nDOF,Mesh.nDOF);

%% insert the edge contributions

%% the edge interpolation matrix, from nodes to Gauss points
lambda = sqrt(15)/10;  %% the sign of lambda is flipped, compared to notes
MB = [1 lambda lambda^2 lambda^3;1 0 0 0;1 -lambda lambda^2 -lambda^3]...
     *[-1 9 9 -1;2 -54 +54 -2;+36 -36 -36 +36;-72 +216 -216 +72]/16;

for k = 1:size(Mesh.edges,1)
  if Mesh.edgesT(k)<-1  % it is a Neumann edge
    cor = Mesh.nodes(Mesh.edges(k,:),:); % the four nodes on the edge
    p2 =  (cor(1,:)+cor(4,:))/2; % the three Gauss points on the edge
    p3 =  cor(4,:); 
    L = norm(p3-p2)/9;  % length of edge, divided by 18
    vec_diff = sqrt(0.6)*(p3-p2);
    p1 = p2 + vec_diff;  %% the Gauss points on the edge
    p3 = p2 - vec_diff;
    
    if gN1Func == 0
      edgeVec = zeros(4,1);
    else
      if ischar(gN1Func)
	gN1 = feval(gN1Func,[p1;p2;p3]);
	edgeVec = L*MB'*(gN1.*[5;8;5]);
      else
	edgeVec = L*MB'*(gN1Func*[5;8;5]);
      endif
    endif% gN1Func

    if gN2Func == 0
      B = zeros(4,4);
    else
      if ischar(gN2Func)
	gN2 = feval(gN2Func,[p1;p2;p3]);
	B = L*MB'*diag(gN2.*[5;8;5])*MB;
      else
	B = L*MB'*diag(gN2Func.*[5;8;5])*MB;
      endif
    endif% gN2func

    if gDFunc == 0  %% evaluate Dirichlet values
       gD = zeros(4,1);
    elseif ischar(gDFunc) gD = feval(gDFunc,cor);
    else                  gD = gDFunc*ones(4,1);
    endif% gDFunc

    dofs = Mesh.node2DOF(Mesh.edges(k,:)); %% mid points 2&3 are certainly free
    if dofs(1)>0   %% node 1 free
      if dofs(4)>0 %% both end points free
	gVec(dofs)      -= edgeVec;
        gMat(dofs,dofs) -= B;
      else  % node 1 free, node 4 Dirichlet
	gVec(dofs(1:3)) -= edgeVec(1:3) + B(1:3,4)*gD(4);
	gMat(dofs(1:3),dofs(1:3)) -= B(1:3,1:3);
      endif% dofs(4)>0
    elseif dofs(4)>0  %% node 1 Dirichlet, node 4 free
      gVec(dofs(2:4)) -= edgeVec(2:4) + B(2:4,1)*gD(1);
      gMat(dofs(2:4),dofs(2:4)) -= B(2:4,2:4);
    endif % dofs(1)>0
  endif % Neumann edge
endfor  % k edges
endfunction
