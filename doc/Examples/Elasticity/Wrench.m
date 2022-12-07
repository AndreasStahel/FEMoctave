load WrenchData.m                 %% load the contour data
scale = 0.15/max(x);              %% scale the contour data
x = scale*x; y = scale*y;

Order = 3;                       %% select the order of the elements 1,2 or 3
BC = -22*ones(size(x));          %% default is a force free boundary
BC([1 4]) = -11; BC(17) = -23;   %% fixed at the two horizontal section on the left
                                 %% vertical force on top right segment
Load = 100/(0.05*0.005);         %% 100 N, distributed over length 0.05 and width 0.005

%% create the mesh
Mesh = CreateMeshTriangle('Wrench',[x,y,BC],0.01^2/4);
switch Order
  case 2    Mesh = MeshUpgrade(Mesh,'quadratic');
  case 3    Mesh = MeshUpgrade(Mesh,'cubic');
endswitch

E = 200e9; nu = 0.25; gN = {0,-Load};              %% data for steel
[u1,u2] = PlaneStress(Mesh,E,nu,{0,0},{0,0},gN);   %% solve the plane stress problem

%% display the original and deformed wrench, with the applied force
scale = 0.001*max(y)/max(u2);
x_force = linspace(x(17),x(18),8); y_force = 0.038*ones(size(x_force))+0.02;
vec_x = zeros(size(x_force)); vec_y = -0.02*ones(size(x_force));
figure(1); clf
trimesh(Mesh.elem,Mesh.nodes(:,1),Mesh.nodes(:,2),'color','green','linewidth',1)
hold on
trimesh(Mesh.elem,Mesh.nodes(:,1)+scale*u1,Mesh.nodes(:,2)+scale*u2,'color','red','linewidth',1)
quiver(x_force,y_force,vec_x,vec_y,0)
hold off; axis equal; xlim([-0.01, 0.16])

[sigma_x,sigma_y,tau_xy] = EvaluateStress(Mesh,u1,u2,E,nu);  %% basic stress
vonMises = EvaluateVonMises(sigma_x,sigma_y,tau_xy);         %% von Mises stress

xi = linspace(0.05,0.15,101)'; yi = interp1(x(14:19),y(14:19),xi);
sigma_y_interp = FEMgriddata(Mesh,sigma_y,xi,yi);

figure(2); plot(xi,sigma_y_interp/1e6)
           xlabel('x'); ylabel('\sigma_y [MPa]'); xlim([0.05,0.15])

figure(3); clf; FEMtrimesh(Mesh,vonMises/1e6)
            zlabel('von Mises stress [MPa]');
            colorbar(); view([40 75])
            xlim([0 0.15]); ylim([-0.025 0.09])
            set(gca, 'XTickLabel', [], 'yTickLabel', [], 'zTickLabel', [])

MaxVonMises = max(vonMises); MinVonMises = min(vonMises);
Max_Min_vonMises_MPa = [MaxVonMises,MinVonMises]/1e6
MaxInd = find(vonMises == MaxVonMises); MaxPosition = Mesh.nodes(MaxInd,:);
MinInd = find(vonMises == MinVonMises); MinPosition = Mesh.nodes(MinInd,:);

figure(4); clf; FEMtricontour(Mesh,vonMises/1e6,41)
            hold on; plot([x;x(1)],[y;y(1)],'k');
            plot(MaxPosition(1),MaxPosition(2),'*r',MinPosition(1),MinPosition(2),'*b');
            hold off; axis equal



