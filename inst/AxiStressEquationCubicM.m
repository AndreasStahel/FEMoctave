function [gMat,gVec] = AxiStressEquationCubicM(Mesh,EFunc,nuFunc,fFunc,gDFunc,gNFunc)
%%  [gMat,gVec] = AxiStressEquationCubicM(Mesh,EFunc,nuFunc,fFunc,gDFunc,gNFunc)
%%
%%  setup the equation for an axisymmetric problem with cubic elements

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

w1 = (155 - sqrt(15))/2400;  w2 = (155 + sqrt(15))/2400;  w3 = 0.1125;
w = 2*[w1,w1,w1,w2,w2,w2,w3]'; DW = diag(w);
rV = reshape(Mesh.GP(:,1),nGP/nElem,nElem);  %% radii at Gauss points
den = 1./((1+nuV).*(1-2*nuV));  %% a common denominator
a1 = DW*(rV.*EV.*(1-nuV).*den);
a2 = DW*(EV.*(1-nuV)./rV.*den);
a3 = DW*(rV.*EV.*nuV.*den);
a4 = DW*(EV.*nuV.*den);
a5 = DW*(rV.*EV./(2*(1+nuV)));

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
Si  = zeros(4*100*nElem,1); Sj = Si; Sval = Si; %% maximal number of contributions
gVec = zeros(sum(nDOF),1);

l1 = (12 - 2*sqrt(15))/21;   l2 = (12 + 2*sqrt(15))/21;

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
  Gr = [(cor(3,2)-cor(1,2))*Mxi + (cor(1,2)-cor(2,2))*Mnu]/(2*area);
  Gz = [(cor(1,1)-cor(3,1))*Mxi + (cor(2,1)-cor(1,1))*Mnu]/(2*area);
  mat1 = area*[Gr'*diag(a1(:,k))*Gr+M'*diag(a2(:,k))*M + Gr'*diag(a4(:,k))*M...
       + M'*diag(a4(:,k))*Gr + Gz'*diag(a5(:,k))*Gz,...
      Gr'*diag(a3(:,k))*Gz + M'*diag(a4(:,k))*Gz + Gz'*diag(a5(:,k))*Gr];
  mat2 = area*...
	 [Gz'*diag(a3(:,k))*Gr + Gz'*diag(a4(:,k))*M + Gr'*diag(a5(:,k))*Gz,...
	  Gz'*diag(a1(:,k))*Gz + Gr'*diag(a5(:,k))*Gr];  
  vec1 = area*M'*(w.*f1V(:,k));   vec2 = area*M'*(w.*f2V(:,k));
  dofs1 = Mesh.node2DOF(Mesh.elem(k,:),1);
  dofs2 = Mesh.node2DOF(Mesh.elem(k,:),2);
  for k1 = 1:10
    %%%%%%%%%%% k1 is free for u1  %%%%%%%%%%%%%%%%
    if dofs1(k1)>0 % k1 is free node for u1
      gVec(dofs1(k1))      -= vec1(k1);
      for k2 = 1:10
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
	  Sval(ptrDOF) = mat1(k1,k2+10);  ptrDOF++;
	else  %% k2 is a Dirichlet node for u2
	  if   ischar(gDFunc{2}) gD2 = feval(gDFunc{2},cor(k2,:));
	  else                   gD2 = gDFunc{2};
	  endif% ischchar
 	  gVec(     dofs1(k1)) += mat1(k1,k2+10)*gD2;
	endif %%dofs1(k2)
      endfor  %% k2 = 1:10
    endif %% dofs1(k1)>0

    %%%%%%%%%%% k1 is free for u2  %%%%%%%%%%%%%%%%
    if dofs2(k1)>0 % k1 is free node for u2
      gVec(nDOF(1)+dofs2(k1))      -= vec2(k1);
      for k2 = 1:10
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
	  Sval(ptrDOF) = mat2(k1,k2+10);  ptrDOF++;
	else  %% k2 is a Dirichlet node for u2
	  if   ischar(gDFunc{2}) gD2 = feval(gDFunc{2},cor(k2,:));
	  else                   gD2 = gDFunc{2};
	  endif% ischchar
 	  gVec(nDOF(1)+dofs2(k1)) += mat2(k1,k2+10)*gD2;
	endif %%dofs1(k2)
      endfor %% k2 = 1:10
    endif %% dofs2(k1)>0
  endfor %% k1 = 1:10

endfor % k (elements)

%% add up to create the sparse matrix
%%Sval(find(abs(Sval)<1e-15)) = 0;  % eliminate very small entries
Si = Si(1:ptrDOF-1); Sj = Sj(1:ptrDOF-1); Sval = Sval(1:ptrDOF-1);
gMat = sparse(Si,Sj,Sval,sum(nDOF),sum(nDOF));

%% insert the edge contributions
lambda = sqrt(15)/10;
Mb = [1 -lambda lambda^2 -lambda^3;1 0 0 0;1 lambda lambda^2 lambda^3]...
     *[-1 9 9 -1;2 -54 +54 -2;+36 -36 -36 +36;-72 +216 -216 +72]/16;
Mbc = 1/18*Mb'*diag([5 8 5]);  %% I have no idea why the /2 ????
                               %% TODO: a Comsol computation confirms the result
for k = 1:size(Mesh.edges,1)
  EdgeType = Mesh.edgesT(k);
  EdgeType_x = fix(EdgeType/10);  EdgeType_y = mod(EdgeType,-10);
  if ((EdgeType_x==-3)||(EdgeType_y==-3))
    cor = Mesh.nodes(Mesh.edges(k,:),:); % the three nodes on the edge
    p2 =  (cor(1,:)+cor(4,:))/2; % the three Gauss points on the edge
    p3 =  cor(4,:); 
    vec_diff = sqrt(0.6)*(p3-p2);
    p1 = p2 - vec_diff;  %% the Gauss points on the edge
    p3 = p2 + vec_diff;
    L = norm(cor(4,:)-cor(1,:));
    dofs = Mesh.node2DOF(Mesh.edges(k,:),:);
    dofs1 = dofs(:,1);  dofs2 = dofs(:,2);
    if EdgeType_x == -3  % nonzero force in x-direction
      if ischar(gNFunc{1}) g1 =feval(gNFunc{1},[p1;p2;p3]).*[p1(1);p2(1);p3(1)];
      else                 g1 = gNFunc{1}*[p1(1);p2(1);p3(1)];
      endif
      edgeVec1 = L*Mbc*g1;
      if (dofs1(1)>0)&&(dofs1(4)>0)        %% both end points free
	gVec(dofs1)           -= edgeVec1;
      elseif (dofs1(1)>0)&&(dofs1(4)==0)   %% node 1 free, node 4 Dirichlet
	gVec(dofs1([1,2,3]))  -= edgeVec1([1,2,3]);
      elseif (dofs1(1)==0)&&(dofs1(4)>0)   %% node 1 Dirichlet, node 4 free
	gVec(dofs1([2,3,4]))  -= edgeVec1([2,3,4]);
      endif    
    endif  % EdgeType_x

    if EdgeType_y == -3  % nonzero force in y-direction
      if ischar(gNFunc{2}) g2 =feval(gNFunc{2},[p1;p2;p3]).*[p1(1);p2(1);p3(1)];
      else                 g2 = gNFunc{2}*[p1(1);p2(1);p3(1)];
      endif
      edgeVec2 = L*Mbc*g2;
      if (dofs2(1)>0)&&(dofs2(4)>0)        %% both end points free
	gVec(nDOF(1)+dofs2)           -= edgeVec2;
      elseif (dofs2(1)>0)&&(dofs2(4)==0)   %% node 1 free, node 4 Dirichlet
	gVec(nDOF(1)+dofs2([1,2,3]))  -= edgeVec2([1,2,3]);
      elseif (dofs2(1)==0)&&(dofs2(4)>0)   %% node 1 Dirichlet, node 4 free
	gVec(nDOF(1)+dofs2([2,3,4]))  -= edgeVec2([2,3,4]);
      endif
    endif  % EdgeType_y
  endif %% EdgeType_x||EdgeType_y
endfor  % k edges
endfunction
