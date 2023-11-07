function [gMat,gVec] = AxiStressEquationM(Mesh,EFunc,nuFunc,fFunc,gDFunc,gNFunc)
%%  [gMat,gVec] = AxiStressEquationM(Mesh,EFunc,nuFunc,fFunc,gDFunc,gNFunc)
%%
%%  setup the equation for an axisymmetric problem with linear elements

%% evaluate the functions
nElem = size(Mesh.elem,1); nGP = size(Mesh.GP,1);
nDOF = Mesh.nDOF;

if ischar(EFunc)
  EV = reshape(feval(EFunc,Mesh.GP,Mesh.GPT),nGP/nElem,nElem);
elseif isscalar(EFunc)
  EV = EFunc*ones(nGP/nElem,nElem);
else
  EV = reshape(EFunc,nGP/nElem,nElem);
endif

if ischar(nuFunc)
  nuV = reshape(feval(nuFunc,Mesh.GP,Mesh.GPT),nGP/nElem,nElem);
elseif isscalar(nuFunc)
  nuV = nuFunc*ones(nGP/nElem,nElem);
else
  nuV = reshape(nuFunc,nGP/nElem,nElem);
endif

rV = reshape(Mesh.GP(:,1),nGP/nElem,nElem);  %% radii at Gauss points
den = 1./((1+nuV).*(1-2*nuV));  %% a common denominator
a1 = sum(rV.*EV.*(1-nuV).*den)/3;
a2 = EV.*(1-nuV)./rV.*den/3;
a3 = sum(rV.*EV.*nuV.*den)/3;
a4 = EV.*nuV.*den/3;
a5 = sum(rV.*EV./(2*(1+nuV)))/3;

if ischar(fFunc{1})
  f1V = reshape(feval(fFunc{1},Mesh.GP,Mesh.GPT),nGP/nElem,nElem);
elseif isscalar(fFunc{1})
  f1V = fFunc{1}*ones(nGP/nElem,nElem);
else
  f1V = reshape(fFunc{1},nGP/nElem,nElem);
endif
f1V = f1V.*rV;  %% multiply the values of f with the radii

if ischar(fFunc{2})
  f2V = reshape(feval(fFunc{2},Mesh.GP,Mesh.GPT),nGP/nElem,nElem);
elseif isscalar(fFunc{2})
  f2V = fFunc{2}*ones(nGP/nElem,nElem);
else
  f2V = reshape(fFunc{2},nGP/nElem,nElem);
endif
f2V = f2V.*rV;  %% multiply the values of f with the radii

%% create memory for the sparse matrix and the RHS vector
Si   = zeros(4*9*nElem,1); Sj = Sval = Si;  %% maximal number of contributions
gVec = zeros(sum(nDOF),1);

M = [4 1 1; 1 4 1; 1 1 4]/6;  %% interpolation matrix, is symmetric
% insert the element matrices and vectors into the global matrix
ptrDOF = 1;  %% counter for the DOF we are working on
for k = 1:nElem   %%for each element
  cor = Mesh.nodes(Mesh.elem(k,:),:);  % coordinates of the nodes
 %% compute element stiffness matrix and vector
  area = Mesh.elemArea(k);  % area = 0.5*det(T)
  Gr = -[cor(3,2)-cor(2,2),cor(1,2)-cor(3,2),cor(2,2)-cor(1,2)]/(2*area);
  Gz = -[cor(2,1)-cor(3,1),cor(3,1)-cor(1,1),cor(1,1)-cor(2,1)]/(2*area);
  mat1 = area*[a1(k)*Gr'*Gr + M*diag(a2(:,k))*M + (a4(:,k)*Gr)'*M ...
	       + M*a4(:,k)*Gr + a5(k)*Gz'*Gz ,...
               a3(k)*Gr'*Gz + M*a4(:,k)*Gz + a5(k)*Gz'*Gr];
  mat2 = area*[a3(k)*Gz'*Gr + (a4(:,k)*Gz)'*M + a5(k)*Gr'*Gz,...
	       a1(k)*Gz'*Gz + a5(k)*Gr'*Gr];
  vec1 = area/3*M*f1V(:,k);  vec2 = area/3*M*f2V(:,k);
  dofs1 = Mesh.node2DOF(Mesh.elem(k,:),1);
  dofs2 = Mesh.node2DOF(Mesh.elem(k,:),2);
  for k1 = 1:3
    %%%%%%%%%%% k1 is free for u1  %%%%%%%%%%%%%%%%
    if dofs1(k1)>0 % k1 is free node for u1
      gVec(dofs1(k1))      -= vec1(k1);
      for k2 = 1:3
 	%% gMat(dofs(k1),dofs(k2)) = gMat(dofs(k1),dofs(k2)) + mat(k1,k2);
	if dofs1(k2)>0  % k2 is free node for u1
	  Si(ptrDOF)   = dofs1(k1);   Sj(ptrDOF) = dofs1(k2);
	  Sval(ptrDOF) = mat1(k1,k2);  ptrDOF++;
	else  %% k2 is a Dirichlet node for u1
	  if   ischar(gDFunc{1}) gD1 = feval(gDFunc{1},cor(k2,:));
	  else                   gD1 = gDFunc{1};
	  endif% ischchar
 	  gVec(dofs1(k1)) += mat1(k1,k2)*gD1;
	endif %%dofs1(k2)

	if dofs2(k2)>0  % k2 is free node for u2
	  Si(ptrDOF)   = dofs1(k1); Sj(ptrDOF) = nDOF(1)+dofs2(k2);
	  Sval(ptrDOF) = mat1(k1,k2+3);  ptrDOF++;
	else  %% k2 is a Dirichlet node for u2
	  if   ischar(gDFunc{2}) gD2 = feval(gDFunc{2},cor(k2,:));
	  else                   gD2 = gDFunc{2};
	  endif% ischchar
 	  gVec(     dofs1(k1)) += mat1(k1,k2+3)*gD2;
	endif %%dofs1(k2)
      endfor  %% k2 = 1:3
    endif %% dofs1(k1)>0

    %%%%%%%%%%% k1 is free for u2  %%%%%%%%%%%%%%%%
    if dofs2(k1)>0 % k1 is free node for u2
      gVec(nDOF(1)+dofs2(k1))      -= vec2(k1);
      for k2 = 1:3
 	%% gMat(dofs(k1),dofs(k2)) = gMat(dofs(k1),dofs(k2)) + mat(k1,k2);
	if dofs1(k2)>0  % k2 is free node for u1
	  Si(ptrDOF)   = nDOF(1)+dofs2(k1); Sj(ptrDOF) = dofs1(k2);
	  Sval(ptrDOF) = mat2(k1,k2);  ptrDOF++;
	else  %% k2 is a Dirichlet node for u1
	  if   ischar(gDFunc{1}) gD1 = feval(gDFunc{1},cor(k2,:));
	  else                   gD1 = gDFunc{1};
	  endif% ischchar
 	  gVec(nDOF(1)+dofs2(k1)) += mat2(k1,k2)*gD1;
	endif %%dofs1(k2)

	if dofs2(k2)>0  % k2 is free node for u2
	  Si(ptrDOF)   = nDOF(1)+dofs2(k1); Sj(ptrDOF) = nDOF(1)+dofs2(k2);
	  Sval(ptrDOF) = mat2(k1,k2+3);  ptrDOF++;
	else  %% k2 is a Dirichlet node for u2
	  if   ischar(gDFunc{2}) gD2 = feval(gDFunc{2},cor(k2,:));
	  else                   gD2 = gDFunc{2};
	  endif% ischchar
 	  gVec(nDOF(1)+dofs2(k1)) += mat2(k1,k2+3)*gD2;
	endif %%dofs1(k2)
      endfor %% k2 = 1:3
    endif %% dofs2(k1)>0
  endfor %% k1 = 1:3

endfor % k (elements)

%% add up to create the sparse matrix
%% Sval(find(abs(Sval)<1e-15)) = 0;  % eliminate very small entries
Si = Si(1:ptrDOF-1); Sj = Sj(1:ptrDOF-1); Sval = Sval(1:ptrDOF-1);
gMat = sparse(Si,Sj,Sval,sum(nDOF),sum(nDOF));

alpha = (1-1/sqrt(3))/2;  %% for interpolation along edge
edgeMat =  [(1-alpha), alpha;  alpha , (1-alpha)];

%% insert the edge contributions
for k = 1:size(Mesh.edges,1)
  EdgeType = Mesh.edgesT(k);
  EdgeType_x = fix(EdgeType/10);  EdgeType_y = mod(EdgeType,-10);
  if ((EdgeType_x==-3)||(EdgeType_y==-3))
    cor = Mesh.nodes(Mesh.edges(k,:),:);
    %%    edgeVec = ElementContributionEdge(corners,gNFunc);
    p1 = (cor(1,:)+cor(2,:))/2 - (cor(2,:)-cor(1,:))/(2*sqrt(3));
    p2 = (cor(1,:)+cor(2,:))/2 + (cor(2,:)-cor(1,:))/(2*sqrt(3));
    L = norm(cor(2,:)-cor(1,:))/2;  %% use half of the length
    dofs = Mesh.node2DOF(Mesh.edges(k,:),:);
    dofs1 = dofs(:,1);  dofs2 = dofs(:,2);
    if EdgeType_x == -3  % nonzero force in x-direction
      if ischar(gNFunc{1}) g1 = feval(gNFunc{1},[p1;p2]).*[p1(1);p2(1)];
      else                 g1 = gNFunc{1}*[p1(1);p2(1)];
      endif
      edgeVec1 = L*edgeMat*g1;
      if (dofs1(1)>0)&&(dofs1(2)>0) %% both points free
	gVec(dofs1)         -= edgeVec1;
      elseif (dofs1(1)>0)&&(dofs1(2)==0)   %% node 1 free, node 2 Dirichlet
	gVec(dofs1(1))      -= edgeVec1(1);
      elseif (dofs1(1)==0)&&(dofs1(2)>0)   %% node 1 Dirichlet, node 2 free
	gVec(dofs1(2))      -= edgeVec1(2);
      endif
    endif  % EdgeType_x
    if EdgeType_y == -3  % nonzero force in y-direction
      if ischar(gNFunc{2}) g2 = feval(gNFunc{2},[p1;p2]).*[p1(1);p2(1)];
      else                 g2 = gNFunc{2}.*[p1(1);p2(1)];
      endif
      edgeVec2 = L*edgeMat*g2;
      if (dofs2(1)>0)&&(dofs2(2)>0) %% both points free
	gVec(nDOF(1)+dofs2)    -= edgeVec2;
      elseif (dofs2(1)>0)&&(dofs2(2)==0)   %% node 1 free, node 2 Dirichlet
	gVec(nDOF(1)+dofs2(1)) -= edgeVec2(1);
      elseif (dofs2(1)==0)&&(dofs2(2)>0)   %% node 1 Dirichlet, node 2 free
	gVec(nDOF(1)+dofs2(2)) -= edgeVec2(2);
      endif
    endif  % EdgeType_y
  endif % Neumann edge with nonzero force
endfor  % k edges

endfunction
