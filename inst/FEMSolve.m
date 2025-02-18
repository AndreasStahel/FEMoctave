function u = FEMSolve(FEMmesh,A,b,gD)
%[...] = FEMSolve (...)
%  solves the system of linear equations for a numerical solution of a PDE
%
%  u = FEMSolve(FEMmesh,A,b,gD)
%
% nodes contains information about the mesh
%         see ReadMesh() for the description of the format
% A    is the matrix of the system to be solved.
%      It is stored as a full matrix
% b    is the RHS of the system to be solved.
% 'gD' is the function describing the Dirichlet boundary condition
%
% u    is the vector with the values of the solution
%
%see also ReadMesh, ShowMesh, FEMEquation, ShowSolution, FEMValue

if (nargin!=4) 
%%  help("FEMSolve");
  print_usage();
endif

%% solve the system of linear equations
ug = -A\b;
n  = size(FEMmesh.node2DOF,1);
u  = zeros(n,1);

%% evaluate the function on the Dirichlet section of the boundary and
%% create a vector u with the solution
%%for k = 1:n
%%  dof = FEMmesh.node2DOF(k);
%%  if dof>0
%%    u(k) = ug(dof);
%%  else
%%    if ischar(gD)
%%      u(k) = feval(gD,FEMmesh.nodes(k,:));
%%    else
%%      u(k) = gD;
%%    endif % ischar
%%  endif% dof>0
%%endfor
%% without loops
ind_free = find(FEMmesh.node2DOF>0);
u(ind_free) = ug;
ind_Dirichlet = find(FEMmesh.node2DOF==0);
if ischar(gD)
  u(ind_Dirichlet) = feval(gD,FEMmesh.nodes(ind_Dirichlet,:));
else
  u(ind_Dirichlet) = gD;
endif % ischar
endfunction
