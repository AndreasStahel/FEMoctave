function [gMat,wMat,OneMat] = PStressEquationWM(Mesh,EFunc,nuFunc,wFunc)
%%  [gMat,wMat,OneMat] = PStressEquationWM(Mesh,EFunc,nuFunc,wFunc)
%%
%%  setup the equation for a plane stress problem with linear elements

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

a1 = sum(EV./(1-nuV.^2))/3;
a2 = sum(nuV.*EV./(1-nuV.^2))/3;
a3 = sum(EV./(1+nuV))/6;

if ischar(wFunc)
  wV = reshape(feval(wFunc,Mesh.GP,Mesh.GPT),nGP/nElem,nElem);
elseif isscalar(wFunc)
  wV = wFunc*ones(nGP/nElem,nElem);
else
  wV = reshape(wFunc,nGP/nElem,nElem);
endif

%% create memory for the sparse matrix and the RHS vector
Si   = zeros(4*9*nElem,1); Sj = Sval = Si;  %% maximal number of contributions
Wi = Wj = Wval = Si;
Onei = Onej = Oneval = Si;

M = [4 1 1; 1 4 1; 1 1 4]/6;  %% interpolation matrix, is symmetric
% insert the element matrices and vectors into the global matrix
ptrDOF = 1;  %% counter for the DOF we are working on
ptrW   = 1;  %% counter for the weight matrix
for k = 1:nElem   %%for each element
  cor = Mesh.nodes(Mesh.elem(k,:),:);  % coordinates of the nodes
 %% compute element stiffness matrix and vector
  area = Mesh.elemArea(k);  % area = 0.5*det(T)
  Gx = -[cor(3,2)-cor(2,2),cor(1,2)-cor(3,2),cor(2,2)-cor(1,2)]/(2*area);
  Gy = -[cor(2,1)-cor(3,1),cor(3,1)-cor(1,1),cor(1,1)-cor(2,1)]/(2*area);
  mat1 = area*[a1(k)*Gx'*Gx + a3(k)*Gy'*Gy , a2(k)*Gx'*Gy + a3(k)*Gy'*Gx];
  mat2 = area*[a2(k)*Gy'*Gx + a3(k)*Gx'*Gy , a1(k)*Gy'*Gy + a3(k)*Gx'*Gx];
  wMat = area/3*M*diag(wV(:,k))*M;
  OneMat = area/3*M*M;
  dofs1 = Mesh.node2DOF(Mesh.elem(k,:),1);
  dofs2 = Mesh.node2DOF(Mesh.elem(k,:),2);
  for k1 = 1:3
    %%%%%%%%%%% k1 is free for u1  %%%%%%%%%%%%%%%%
    if dofs1(k1)>0 % k1 is free node for u1
      for k2 = 1:3
 	%% gMat(dofs(k1),dofs(k2)) = gMat(dofs(k1),dofs(k2)) + mat(k1,k2);
	if dofs1(k2)>0  % k2 is free node for u1
	  Wi(ptrW)     = dofs1(k1);   Wj(ptrW) = dofs1(k2);
	  Wval(ptrW)   = wMat(k1,k2);
	  Onei(ptrW)   = dofs1(k1);   Onej(ptrW) = dofs1(k2);
	  Oneval(ptrW) = OneMat(k1,k2);  ptrW++;
	  Si(ptrDOF)   = dofs1(k1);   Sj(ptrDOF) = dofs1(k2);
	  Sval(ptrDOF) = mat1(k1,k2);  ptrDOF++;
	endif %%dofs1(k2)
	
	if dofs2(k2)>0  % k2 is free node for u2
	  Si(ptrDOF)   = dofs1(k1); Sj(ptrDOF) = nDOF(1)+dofs2(k2);
	  Sval(ptrDOF) = mat1(k1,k2+3);  ptrDOF++;
	endif %%dofs1(k2)
      endfor  %% k2 = 1:3
    endif %% dofs1(k1)>0

    %%%%%%%%%%% k1 is free for u2  %%%%%%%%%%%%%%%%
    if dofs2(k1)>0 % k1 is free node for u2
      for k2 = 1:3
 	%% gMat(dofs(k1),dofs(k2)) = gMat(dofs(k1),dofs(k2)) + mat(k1,k2);
	if dofs1(k2)>0  % k2 is free node for u1
	  Si(ptrDOF)   = nDOF(1) + dofs2(k1); Sj(ptrDOF) = dofs1(k2);
	  Sval(ptrDOF) = mat2(k1,k2);  ptrDOF++;
	endif %%dofs1(k2)
	
	if dofs2(k2)>0  % k2 is free node for u2
	  Wi(ptrW)     = nDOF(1)+dofs2(k1); Wj(ptrW) = nDOF(1)+dofs2(k2);
	  Wval(ptrW)   = wMat(k1,k2);
	  Onei(ptrW)   = nDOF(1)+dofs2(k1); Onej(ptrW) = nDOF(1)+dofs2(k2);
	  Oneval(ptrW) = OneMat(k1,k2);  ptrW++;
	  Si(ptrDOF)   = nDOF(1)+dofs2(k1); Sj(ptrDOF) = nDOF(1)+dofs2(k2);
	  Sval(ptrDOF) = mat2(k1,k2+3);  ptrDOF++;
	endif %%dofs1(k2)
      endfor %% k2 = 1:3
    endif %% dofs2(k1)>0
  endfor %% k1 = 1:3

endfor % k (elements)

%% add up to create the sparse matrix
%% Sval(find(abs(Sval)<1e-15)) = 0;  % eliminate very small entries
sum_nDOF = sum(nDOF);
Si = Si(1:ptrDOF-1); Sj = Sj(1:ptrDOF-1); Sval = Sval(1:ptrDOF-1);
gMat = sparse(Si,Sj,Sval,sum_nDOF,sum_nDOF);
Wi = Wi(1:ptrW-1); Wj = Wj(1:ptrW-1); Wval = Wval(1:ptrW-1);
wMat = sparse(Wi,Wj,Wval,sum_nDOF,sum_nDOF);
Onei = Onei(1:ptrW-1); Onej = Onej(1:ptrW-1); Oneval = Oneval(1:ptrW-1);
OneMat = sparse(Onei,Onej,Oneval,sum_nDOF,sum_nDOF);
endfunction
