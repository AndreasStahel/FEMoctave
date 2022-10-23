function [gMat,gVec] = FEMEquationQuadM(Mesh,aFunc,bFunc,fFunc,gDFunc,gN1Func,gN2Func)
%[...] = FEMEquationQuadM (...)
%  sets up the system of linear equations for a numerical solution of a PDE
%
%  [A,b,n2d] = FEMEquationQuadM(Mesh,'a','b','f','gD','gN1','gN2')
%    Mesh is the mesh describing the domain\n\
%         see ReadMesh() for the description of the format
%   'a','b','f','gD','gN1','gN2' are the functions and coefficients 
%         for the boundary value problem. They can be given as a scalar value
%         or as a sting with the function name 
%
%  -div(a*grad u) + b*u = f            in domain
%                     u = gD           on Dirichlet section of the boundary
%               a*du/dn = gN1+gN2*u    on Neumann section of the boundary 
%
%
% A   is the matrix of the system to be solved.
% b   is the RHS of the system to be solved.
%

%% evaluate the functions a b and f

nElem = size(Mesh.elem,1);
nGP   = size(Mesh.GP,1)  ;

if ischar(aFunc)
  aV = reshape(feval(aFunc,Mesh.GP,Mesh.GPT),nGP/nElem,nElem);
elseif isscalar(aFunc)
  aV = aFunc*ones(nGP/nElem,nElem);
else
  aV = reshape(aFunc,nGP/nElem,nElem);
endif

if ischar(bFunc)
  bV = reshape(feval(bFunc,Mesh.GP,Mesh.GPT),nGP/nElem,nElem);
elseif isscalar(bFunc)
  bV = bFunc*ones(nGP/nElem,nElem);
else
  bV = reshape(bFunc,nGP/nElem,nElem);
endif

if ischar(fFunc)
  fV = reshape(feval(fFunc,Mesh.GP,Mesh.GPT),nGP/nElem,nElem);
elseif isscalar(fFunc)
  fV = fFunc*ones(nGP/nElem,nElem);
else
  fV = reshape(fFunc,nGP/nElem,nElem);
endif

%% create memory for the sparse matrix and the RHS vector
Si   = zeros(36*nElem,1); Sj = Si; Sval = Si; %% maximal number of contributions
gVec = zeros(Mesh.nDOF,1);

l1 = (12 - 2*sqrt(15))/21;   l2 = (12 + 2*sqrt(15))/21;
w1 = (155 - sqrt(15))/2400;  w2 = (155 + sqrt(15))/2400;
w3 = 0.1125;         w = [w1,w1,w1,w2,w2,w2,w3]';

xi = [l1/2, 1-l1, l1/2, l2/2, 1-l2, l2/2, 1/3]';
nu = [l1/2, l1/2, 1-l1, l2/2, l2/2, 1-l2, 1/3]';

%% the interpolation matrices for the function and its partial derivatives
M = [(1-xi-nu).*(1-2*xi-2*nu) xi.*(2*xi-1) nu.*(2*nu-1) 4*xi.*nu 4*nu.*(1-xi-nu) 4*xi.*(1-xi-nu)];
Mxi = [-3+4*(xi+nu) 4*xi-1 0*xi 4*nu -4*nu  4-8*xi-4*nu];
Mnu = [-3+4*(xi+nu) 0*xi 4*nu-1 4*xi  4-4*xi-8*nu -4*xi];

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
  mat = 0.5/area*(Ax+Ay) + 2*area*B;
  vec = -2*area*M'*(w.*fV(:,k));
  dofs = Mesh.node2DOF(Mesh.elem(k,:));
  for k1 = 1:6
    if dofs(k1)>0 % k1 is free node
      gVec(dofs(k1)) +=  vec(k1);
      for k2 = 1:6
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
alpha = sqrt(0.6)/2; %% the edge interpolation matrix, transposed
MBtr = [alpha*(+1+2*alpha) 1-4*alpha^2 alpha*(-1+2*alpha);
        0 1 0;
        alpha*(-1+2*alpha) 1-4*alpha^2 alpha*(+1+2*alpha)]';

%MBtr = [alpha*(+1+2*alpha) 0 alpha*(-1+2*alpha);
%        0.5 0 0.5;
%        alpha*(-1+2*alpha) 0 alpha*(+1+2*alpha)]';

for k = 1:size(Mesh.edges)(1)
  if Mesh.edgesT(k)<-1  % it is a Neumann edge
    cor = Mesh.nodes(Mesh.edges(k,:),:); % the three nodes on the edge
    p2 =  cor(2,:); % the three Gauss points on the edge
    p3 =  cor(3,:); 
    L = norm(p3-p2)/9;  % length of edge, divided by 18
    vec_diff = sqrt(0.6)*(p3-p2);
    p1 = p2 - vec_diff;  %% the Gauss points on the edge
    p3 = p2 + vec_diff;

    if gN1Func == 0
      edgeVec = zeros(3,1);
    else
      if ischar(gN1Func)
	g = feval(gN1Func,[p1;p2;p3]);
	edgeVec = L*MBtr*(g.*[5;8;5]);
      else
	edgeVec = L*MBtr*(gN1Func*[5;8;5]);
      endif
    endif% gN1Func

    if gN2Func == 0
      B = zeros(3,3);
    else
      if ischar(gN2Func)
	gN2 = feval(gN2Func,[p1;p2;p3]);
	B   = L*MBtr*diag(gN2.*[5;8;5])*MBtr';
      else
	B = L*MBtr*diag(gN2Func.*[5;8;5])*MBtr';
      endif
    endif% gN2func
    
    if gDFunc == 0  %% evaluate Dirichlet values
      gD = zeros(3,1);
    elseif ischar(gDFunc) gD = feval(gDFunc,cor);
    else                  gD = gDFunc*ones(3,1);
    endif% gDFunc

    dofs = Mesh.node2DOF(Mesh.edges(k,:)); %% mid point is certainly free
    if dofs(1)>0   %% node 1 free
      if dofs(3)>0 %% both end points free
	gVec(dofs) -= edgeVec;
        gMat(dofs,dofs) -= B;
      else  % node 1 free, node 3 Dirichlet
	gVec(dofs(1:2)) -= edgeVec(1:2) + B(1:2,3)*gD(3);
	gMat(dofs(1:2),dofs(1:2)) -= B(1:2,1:2);
      endif% dofs(3)>0
    elseif dofs(3)>0  %% node 1 Dirichlet, node 3 free
      gVec(dofs(2:3)) -= edgeVec(2:3) + B(2:3,1)*gD(1);
      gMat(dofs(2:3),dofs(2:3)) -= B(2:3,2:3);
    endif % dofs(1)>0
  endif % Neumann edge
endfor  % k edges
endfunction
