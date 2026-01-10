## -*- texinfo -*-
## @deftypefn  {} {} HeatSemilinearLoop.m
##
## This is a demo file  inside the `doc/Demos/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

h  = 0.1; %%typical radius of triangles
nodes = [0 0 -1; pi 0, -1; pi pi -1; 0 pi -1];
FEMmesh = CreateMeshTriangle("test",nodes,h^2/2);
FEMmesh = MeshUpgrade(FEMmesh,'cubic');

f = @(xy,t,u) u.*(1-u) + exp(-2*t)*sin(xy(:,1)).^2 .* sin(xy(:,2)).^2;
t_end = 1;
u0 = @(xy)sin(xy(:,1)) .* sin(xy(:,2));
u_exact = exp(-t_end)*sin(FEMmesh.nodes(:,1)) .* sin(FEMmesh.nodes(:,2));

Nsteps = [10,20,40,80,160];
Diff_Max_RMS = zeros(length(Nsteps),4);
tic();
for jj = 1:length(Nsteps)
  Steps = [1,Nsteps(jj)];
  Solver = 'CNPC';
  [u_dyn,t] = IBVP2DNL(FEMmesh,1,1,0,0,0,f,0,0,0,u0,0,t_end,Steps,'solver',Solver);
  u_end = u_dyn(:,end);
  Difference1 = abs(u_end-u_exact);
  Solver = 'CNexp';
  [u_dyn,t] = IBVP2DNL(FEMmesh,1,1,0,0,0,f,0,0,0,u0,0,t_end,Steps,'solver',Solver);
  u_end = u_dyn(:,end);
  Difference2 = abs(u_end-u_exact);
  Diff_Max_RMS(jj,:) = [max(Difference1),sqrt(mean(Difference1.^2)),...
                        max(Difference2),sqrt(mean(Difference2.^2))];
endfor
timing = toc()
h = t_end./Nsteps;
figure(1); plot(h,Diff_Max_RMS(:,1),'+ -',h,Diff_Max_RMS(:,2),'+ -',...
                h,Diff_Max_RMS(:,3),'+ -',h,Diff_Max_RMS(:,4),'+ -')
           xlabel('size of steps'); ylabel('error');
           legend('CNPC max', 'CNPC RMS','CNEXP max','CNEXP RMS','location','northwest');
figure(2); loglog(h,Diff_Max_RMS(:,1),'+ -',h,Diff_Max_RMS(:,2),'+ -',...
                h,Diff_Max_RMS(:,3),'+ -',h,Diff_Max_RMS(:,4),'+ -')
           xlabel('size of steps'); ylabel('error'); xlim([1/200,1/8]);
           legend('CNPC max', 'CNPC RMS','CNEXP max','CNEXP RMS','location','northwest');
