## -*- texinfo -*-
## @deftypefn{function file}{}[@var{u1},@var{u2}] = PlaneStrain(@var{mesh},@var{E},@var{nu},@var{f},@var{gD},@var{gN})
##
##   solve an plane strain problem
##
##@verbatim
##    plane strain equation    in domain
##                   u = gD    on Gamma_1
##       force density = gN    on Gamma_2 
##       force density = 0     on Gamma_3 
##@end verbatim
##
##parameters:
##@itemize
##@item @var{mesh} is the mesh describing the domain and the boundary types
##@item @var{E},@var{nu} Young's modulus and Poisson's ratio for the material
##@item @var{f = @{f1,f2@}} a cell array with the two components of the volume forces
##@item @var{gD = @{gD1,gD2@}} a cell array with the two components of the prescribed displacements on the boundary section Gamma_1
##@item @var{gN = @{gN1,gN2@}} a cell array with the two components of the surface forces on the boundary section Gamma_2
##@item Any constant function can be given by its scalar value
##@item Any function can be given by a string with the function name
##@item The functions @var{E}, @var{nu}, @var{f1} and @var{f2} may also be given as vectors with the values of the function at the Gauss points
##@end itemize
##
##return values
##@itemize
##@item @var{u1}  vector with the values of the x-displacement at the nodes
##@item @var{u2}  vector with the values of the y-displacement at the nodes
##@end itemize
##
## @c Will be cut out in ??? info file and replaced with the same
## @c references explicitly there, since references to core Octave
## @c functions are not automatically transformed from here to there.
## @c BEGIN_CUT_TEXINFO
## @seealso{PlaneStress}
## @c END_CUT_TEXINFO
## @end deftypefn

function [u1,u2] = PlaneStrain(Mesh,E,nu,f,gD,gN)
  if nargin ~=6
    print_usage();
  endif

  nElem = size(Mesh.elem,1); nGP = size(Mesh.GP,1);
  if ischar(E)
  EV = reshape(feval(E,Mesh.GP,Mesh.GPT),nGP/nElem,nElem);
  elseif isscalar(E)
    EV = E*ones(nGP/nElem,nElem);
  else
    EV = reshape(E,nGP/nElem,nElem);
  endif
  
  if ischar(nu)
    nuV = reshape(feval(nu,Mesh.GP,Mesh.GPT),nGP/nElem,nElem);
  elseif isscalar(nu)
    nuV = nu*ones(nGP/nElem,nElem);
  else
    nuV = reshape(nu,nGP/nElem,nElem);
  endif
  Estar = EV./(1-nuV.^2);
  nustar = nuV./(1-nuV);

  switch Mesh.type
    case 'linear' %% first order elements
      [A,b] = PStressEquationM (Mesh,Estar,nustar,f,gD,gN);
    case 'quadratic'  %% second order elements
      [A,b] = PStressEquationQuadM (Mesh,Estar,nustar,f,gD,gN);
    case 'cubic'      %% third order elements
      [A,b] = PStressEquationCubicM (Mesh,Estar,nustar,f,gD,gN);
  endswitch

  
  ug = -A\b;
    nDOF = Mesh.nDOF;   n = size(Mesh.nodesT,1);  u1 = zeros(n,1); u2 = u1;

  ind_free1 = find(Mesh.node2DOF(:,1)>0);
  ind_free2 = find(Mesh.node2DOF(:,2)>0);

  %%u([ind_free1;n+ind_free2]) = ug;
  u1(ind_free1) = ug([1:nDOF(1)]);
  u2(ind_free2) = ug([nDOF(1)+1:end]);
  ind_Dirichlet1 = find(Mesh.node2DOF(:,1)==0);
  ind_Dirichlet2 = find(Mesh.node2DOF(:,2)==0);

  if isscalar(gD{1})
    u1(ind_Dirichlet1)   = gD{1};
  else
    u1(ind_Dirichlet1)   = feval(gD{1},Mesh.nodes(ind_Dirichlet1,:));
  endif
  if isscalar(gD{2})
    u2(ind_Dirichlet2) = gD{2};
  else
    u2(ind_Dirichlet2) = feval(gD{2},Mesh.nodes(ind_Dirichlet2,:));
  endif
endfunction
