%% create a circle with mesh
x_p = [0;1;1;-1;-1;0]; y_p = [ 0;0;1;1;-1;-1];
FEMmesh = CreateMeshTriangle("test",[x_p,y_p,-ones(size(x_p))], 0.01);
FEMmesh = MeshUpgrade(FEMmesh,'quadratic');

function res = gD(xy)
  phi = mod(atan2(xy(:,2),xy(:,1)),2*pi);
  res = (xy(:,1).^2+ xy(:,2).^2).^(1/3).*sin(2/3*phi);
endfunction

u = BVP2Dsym(FEMmesh,1,0,0,'gD',0,0);
figure(1); FEMtrimesh(FEMmesh,u);
           xlabel("x"); ylabel("y"); title('FEM solution')
           view([30,30])

u_exact = gD(FEMmesh.nodes);
figure(2); FEMtrimesh(FEMmesh,-u+u_exact);
           xlabel("x"); ylabel("y"); title('Error of FEM solution')
           view([30,30])


[ux,uy] = FEMEvaluateGradient(FEMmesh,u);
figure(3); FEMtrimesh(FEMmesh,ux);
           xlabel("x"); ylabel("y"); title('FEM solution, u_x')
           view([30,30])

figure(4); FEMtrimesh(FEMmesh,uy);
           xlabel("x"); ylabel("y"); title('FEM solution, u_y')
           view([30,30])


figure(5); FEMtrimesh(FEMmesh,sqrt(ux.^2+uy.^2));
           xlabel("x"); ylabel("y"); title('FEM solution, norm of gradient')
           view([30,30])

