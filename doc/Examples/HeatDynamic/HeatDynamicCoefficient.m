%% parameters
h = 1.2;  l = 0.5; Nt = 60;  %% number of time steps
FEMmesh = CreateMeshTriangle('Test',...
    [0 0,-2; 2+l 0 -2; 2+l 1, -2; 1+l 1 -2; 1+l h -2; 1 h -2; 1 1 -2; 0 1 -1],0.01);
FEMmesh = MeshUpgrade(FEMmesh,'quadratic');
%%FEMmesh = MeshUpgrade(FEMmesh,'cubic');

figure(1); FEMtrimesh(FEMmesh);
           axis equal; xlabel('x'); ylabel('y')

function res = a(xy)
  l = 0.5;
  res = ones(size(xy,1),1);
  res(find(abs(xy(:,1)-1-l/2)<l/2)) *= 1/6; 
endfunction

[u t] = IBVP2D(FEMmesh,1,'a',0, 0, 0, 0,1, 0,  0,  0, 0, 10, [Nt,10]);

figure(2); FEMtrimesh(FEMmesh,u(:,end))
           xlabel('x'); ylabel('y'); zlim([0,1]); view([10 30]); caxis([0,1]);
	   text(0.2,0.2,0.2,sprintf('t = %4.2f',t(end))); zlabel('temperature')
	   
	   
figure(3); FEMtrimesh(FEMmesh,u(:,Nt/2+1))
           xlabel('x'); ylabel('y'); zlim([0,1]); view([10 30]); caxis([0,1])
	   text(0.2,0.2,0.2,sprintf('t = %4.2f',t(Nt/2+1))); zlabel('temperature')

figure(4); FEMtrimesh(FEMmesh,u(:,Nt/3+1))
           xlabel('x'); ylabel('y'); zlim([0,1]); view([10 30]); caxis([0,1])
	   text(0.2,0.2,0.2,sprintf('t = %4.2f',t(Nt/5+1))); zlabel('temperature')


x = linspace(0,2+l,51);  u_int = zeros(size(t,2)-1,size(x,2));
for jj = 2:size(t,2)
  u_int(jj-1,:) = FEMgriddata(FEMmesh,u(:,jj),x,zeros(size(x)));
endfor

figure(10); mesh(x,t(2:end),u_int)
            xlabel('x'); ylabel('t'); zlabel('temperature at y=0')

figure(11); [c,h] = contour(x,t(2:end),u_int,[0:0.1:1]);
            clabel(c,h);
            xlabel('x'); ylabel('t');
	    
figure(12); plot(t(2:end),u_int(:,end))
            xlabel('t'); ylabel('temperature at x=end and y=0')
