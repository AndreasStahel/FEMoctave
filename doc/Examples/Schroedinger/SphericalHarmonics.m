Nphi = 50; Ntheta = 50;
phi = linspace(0,2*pi,Nphi)'; theta = linspace(0,pi,Ntheta)';
Mesh = CreateMeshRect(phi,theta,-2,-2,-1,-1);
Mesh = MeshUpgrade(Mesh,'cubic');

function res = a(PhiTheta)
  Theta = PhiTheta(:,2);
  res = [1./sin(Theta), sin(Theta), zeros(size(Theta))];
endfunction

function res = w(PhiTheta)
  res = sin(PhiTheta(:,2));
endfunction
Nval = 30;
[Eval,Evec] = BVP2Deig(Mesh,'a',0,'w',0,Nval);

flag = zeros(Nval,1);  %% select the function different from zero at midpoint
for ii = 1:Nval
  flag(ii) = abs(FEMgriddata(Mesh,Evec(:,ii),pi,1))<1e-2;
endfor
Indices = find(flag>0);
Eval = Eval(Indices);
Evec = Evec(:,Indices);
Nval = sum(flag);

n = 6;
u = Evec(:,n);
figure(1); FEMtrimesh(Mesh,u)
           xlabel('\phi') ; ylabel('\theta'); zlabel('u')
           xlim([0,2*pi]); ylim([0,pi])

figure(2); clf; FEMtricontour(Mesh,u)
           xlabel('\phi') ; ylabel('\theta');
           xlim([0,2*pi]); ylim([0,pi])

phi = linspace(0,2*pi,100);
for nn = 1:Nval
  u = FEMgriddata(Mesh,Evec(:,nn),phi,ones(size(phi)));
  for m = 1:5
    Corr = corrcoef(u,sin(m*phi));
    if abs(Corr(1,2))>0.95
      disp(['solution ',num2str(nn),' with m = ',num2str(m),...
      ' and eigenvalue = ',num2str(Eval(nn))])
    endif
  endfor
endfor
