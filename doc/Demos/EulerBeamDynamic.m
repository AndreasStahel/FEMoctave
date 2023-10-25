E = 70e9; nu = 0.33; rho = 2.7e3; L = 0.2; H = 0.01;
f = {0,0}; gD = {0,0}; gN = {0,0};
function res = u0Func(xy)
  z = 1.8751; L = max(xy(:,1));
  C = -(cos(z)+cosh(z))/(sin(z)+sinh(z));
  x = z*xy(:,1)/L;
  res = cos(x)-cosh(x) + C*(sin(x)-sinh(x));
  res = -0.1*res./max(abs(res));
endfunction
u0 = {0,'u0Func'}; v0 = {0,0};
t0 = 0; tend = 0.0056; steps = [100,20];

Mesh = CreateMeshRect(linspace(0,L,31),linspace(-H/2,+H/2,3),-22,-22,-11,-22);
Mesh = MeshUpgrade(Mesh,'quadratic');  solver = 'implicit';

[u1_all,u2_all,t] = PlaneStressDynamic(Mesh,E,nu,rho,f,gD,gN,u0,v0,t0,tend,...
                                       steps,'solver',solver);
u1 = u1_all(:,end); u2 = u2_all(:,end);

figure(1); FEMtrimesh(Mesh,u2); xlabel('x'); ylabel('y'); zlabel('u_2');
           ylim([-H/2,H/2])

ind = find((Mesh.nodes(:,1)==L).*(Mesh.nodes(:,2)==0));
u2_t = u2_all(ind,:);
figure(2); plot(t*1e3,u2_t); xlabel('t [ms]'); ylabel('u_2(L,0,t)');
           xlim([0,max(t)*1e3]); ylim(0.11*[-1,+1])

