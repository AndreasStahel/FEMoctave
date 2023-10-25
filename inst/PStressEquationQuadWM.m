function [gMat,wMat,OneMat] = PStressEquationQuadWM(Mesh,EFunc,nuFunc,wFunc)
%%  [gMat,wMat,OneMat] = PStressEquationQuadWM(Mesh,EFunc,nuFunc,wFunc)
%%
%%  setup the equation for a plane stress problem with quadratic elements

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

a1 = EV./(1-nuV.^2);
a2 = nuV.*EV./(1-nuV.^2);
a3 = EV./(1+nuV)/2;

if ischar(wFunc)
  wV = reshape(feval(wFunc,Mesh.GP,Mesh.GPT),nGP/nElem,nElem);
elseif isscalar(wFunc)
  wV = wFunc*ones(nGP/nElem,nElem);
else
  wV = reshape(wFunc,nGP/nElem,nElem);
endif


%% create memory for the sparse matrix and the RHS vector
Si   = zeros(4*36*nElem,1); Sj = Si; Sval = Si; %% maximal number of contributions
Wi = Wj = Wval = Si;
Onei = Onej = Oneval = Si;

l1 = (12 - 2*sqrt(15))/21;   l2 = (12 + 2*sqrt(15))/21;
w1 = (155 - sqrt(15))/2400;  w2 = (155 + sqrt(15))/2400;
w3 = 0.1125;         w = 2*[w1,w1,w1,w2,w2,w2,w3]';

xi = [l1/2, 1-l1, l1/2, l2/2, 1-l2, l2/2, 1/3]';
nu = [l1/2, l1/2, 1-l1, l2/2, l2/2, 1-l2, 1/3]';

%% the interpolation matrices for the function and its partial derivatives
M = [(1-xi-nu).*(1-2*xi-2*nu) xi.*(2*xi-1) nu.*(2*nu-1) 4*xi.*nu 4*nu.*(1-xi-nu) 4*xi.*(1-xi-nu)];
Mxi = [-3+4*(xi+nu) 4*xi-1 0*xi 4*nu -4*nu  4-8*xi-4*nu];
Mnu = [-3+4*(xi+nu) 0*xi 4*nu-1 4*xi  4-4*xi-8*nu -4*xi];

% insert the element matrices and vectors into the global matrix
ptrDOF = 1;    %% counter for the DOF we are working on
ptrW   = 1;    %% counter for the weight matrix
ptrOne = 1;  %% counter for the One matrix
for k = 1:nElem   %%for each element
  cor = Mesh.nodes(Mesh.elem(k,:),:);  % coordinates of the nodes
 %% compute element stiffness matrix and vector
  area = Mesh.elemArea(k);  % area = 0.5*det(T)
  Gx = [(cor(3,2)-cor(1,2))*Mxi + (cor(1,2)-cor(2,2))*Mnu]/(2*area);
  Gy = [(cor(1,1)-cor(3,1))*Mxi + (cor(2,1)-cor(1,1))*Mnu]/(2*area);
  mat1 = [Gx'*diag(w.*a1(:,k))*Gx + Gy'*diag(w.*a3(:,k))*Gy,...
	  Gx'*diag(w.*a2(:,k))*Gy + Gy'*diag(w.*a3(:,k))*Gx]*area;
  mat2 = [Gy'*diag(w.*a2(:,k))*Gx + Gx'*diag(w.*a3(:,k))*Gy,...
	  Gy'*diag(w.*a1(:,k))*Gy + Gx'*diag(w.*a3(:,k))*Gx]*area;
  wMat = area*M'*diag(w.*wV(:,k))*M;
  OneMat = area*M'*diag(w)*M;
  
  dofs1 = Mesh.node2DOF(Mesh.elem(k,:),1);
  dofs2 = Mesh.node2DOF(Mesh.elem(k,:),2);
  for k1 = 1:6
    %%%%%%%%%%% k1 is free for u1  %%%%%%%%%%%%%%%%
    if dofs1(k1)>0 % k1 is free node for u1
      for k2 = 1:6
 	%% gMat(dofs(k1),dofs(k2)) = gMat(dofs(k1),dofs(k2)) + mat(k1,k2);
	if dofs1(k2)>0  % k2 is free node for u1
	  Wi(ptrW)     = dofs1(k1);    Wj(ptrW) = dofs1(k2);
	  Wval(ptrW)   = wMat(k1,k2);
	  Onei(ptrW)   = dofs1(k1);    Onej(ptrW) = dofs1(k2);
	  Oneval(ptrW) = OneMat(k1,k2);  ptrW++;
	  Si(ptrDOF)   = dofs1(k1);   Sj(ptrDOF) = dofs1(k2);
	  Sval(ptrDOF) = mat1(k1,k2);  ptrDOF++;
	endif %%dofs1(k2)
	
	if dofs2(k2)>0  % k2 is free node for u2
	  Si(ptrDOF)   = dofs1(k1); Sj(ptrDOF) = nDOF(1)+dofs2(k2);
	  Sval(ptrDOF) = mat1(k1,k2+6);  ptrDOF++;
	endif %%dofs1(k2)
      endfor  %% k2 = 1:6
    endif %% dofs1(k1)>0

    %%%%%%%%%%% k1 is free for u2  %%%%%%%%%%%%%%%%
    if dofs2(k1)>0 % k1 is free node for u2
      for k2 = 1:6
 	%% gMat(dofs(k1),dofs(k2)) = gMat(dofs(k1),dofs(k2)) + mat(k1,k2);
	if dofs1(k2)>0  % k2 is free node for u1
	  Si(ptrDOF)   = nDOF(1)+dofs2(k1); Sj(ptrDOF) = dofs1(k2);
	  Sval(ptrDOF) = mat2(k1,k2);  ptrDOF++;
	endif %%dofs1(k2)
	
	if dofs2(k2)>0  % k2 is free node for u2
	  Wi(ptrW)     = nDOF(1)+dofs2(k1); Wj(ptrW) = nDOF(1)+dofs2(k2);
	  Wval(ptrW)   = wMat(k1,k2);
	  Onei(ptrW)   = nDOF(1)+dofs2(k1); Onej(ptrW) = nDOF(1)+dofs2(k2);
	  Oneval(ptrW) = OneMat(k1,k2);  ptrW++;
	  Si(ptrDOF)   = nDOF(1)+dofs2(k1); Sj(ptrDOF) = nDOF(1)+dofs2(k2);
	  Sval(ptrDOF) = mat2(k1,k2+6);  ptrDOF++;
	endif %%dofs1(k2)
      endfor %% k2 = 1:6
    endif %% dofs2(k1)>0
  endfor %% k1 = 1:6

endfor % k (elements)

%% add up to create the sparse matrix
%%Sval(find(abs(Sval)<1e-15)) = 0;  % eliminate very small entries
sum_nDOF = sum(nDOF);
Si = Si(1:ptrDOF-1); Sj = Sj(1:ptrDOF-1); Sval = Sval(1:ptrDOF-1);
gMat = sparse(Si,Sj,Sval,sum_nDOF,sum_nDOF);
Wi = Wi(1:ptrW-1); Wj = Wj(1:ptrW-1); Wval = Wval(1:ptrW-1);
wMat = sparse(Wi,Wj,Wval,sum_nDOF,sum_nDOF);
Onei = Onei(1:ptrW-1); Onej = Onej(1:ptrW-1); Oneval = Oneval(1:ptrW-1);
OneMat = sparse(Onei,Onej,Oneval,sum_nDOF,sum_nDOF);
endfunction
