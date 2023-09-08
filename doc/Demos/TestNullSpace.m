%%  script to illustrate the three dimensional null space

Mesh = CreateMeshTriangle('test',[1 0 -22;2 0 -22;2 2 -22; 1 1 -22],0.01);
if 1
  figure(1); FEMtrimesh(Mesh)
endif

E = 1e9; nu = 0.3; f = {0,0}; gD = {0,0}; gN = {0,0};  %% set the parameters
if 0  %% plane stress
  [A,g] = PStressEquationM(Mesh,E,nu,f,gD,gD);           %% determine matrix A
else  %% axially symmetric
  [A,g] = AxiStressEquationM(Mesh,E,nu,f,gD,gD);         %% determine matrix A
endif

A = (A+A')/2;  %% assure that matrix is symmetric, it should be, but rounding errors
EigenValues = eigs(A,6,'sa')           %% find the smallest eigenvalues
%%EigenValues = eigs(A,6,-1e-15)           %% find the smallest eigenvalues

n = size(A,1)/2;
shift_x = [ones(n,1);zeros(n,1)];      %% constant shift in x direction
shift_x = shift_x/norm(shift_x);
Norm_Shift_x  = [norm(A*shift_x),shift_x'*A*shift_x/2]

shift_y = [zeros(n,1);ones(n,1)];      %% constant shift in y direction
shift_y = shift_y/norm(shift_y);
Norm_Shift_y  = [norm(A*shift_y),shift_y'*A*shift_y/2]

%% at point [x,y] add displacement [-y,x]
x = Mesh.nodes(:,1); y = Mesh.nodes(:,2);
rot_vec = [-y;x];                       %% a rotation
rot_vec = rot_vec/norm(rot_vec);
Norm_Rotation = [norm(A*rot_vec),rot_vec'*A*rot_vec/2]

rand_vec = [x.*y;y];    %% an arbitrary displacement vector
rand_vec = rand_vec/norm(rand_vec);
Norm_Random = [norm(A*rand_vec),rand_vec'*A*rand_vec/2]

%% OK, generates the expected result
%% use PlainStressEig()
lambda = PlaneStressEig(Mesh,W,nu,1,6)'
