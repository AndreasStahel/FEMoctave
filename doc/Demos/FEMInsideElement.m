%  illustrate the behavior inside the elements
clear *
MeshType = 'linear'  %% use 'linear', 'quadratic' or 'cubic'
printing = 0;

N = 2;
Mesh = CreateMeshTriangle('test',[0 0 -1;1 0 -1;1 2 -1; 0 1 -1],1/N^2);

switch MeshType
  case 'quadratic'
    Mesh = MeshUpgrade(Mesh,'quadratic');
  case 'cubic'
    Mesh = MeshUpgrade(Mesh,'cubic');
endswitch

figure(1); FEMtrimesh(Mesh)
           xlabel('x'); ylabel('y'); xlim([-0.1,1.1]); ylim([-0.1,2.1])
%axis square
%axis equal
if printing
  print -dpdfcrop doc/InsideElementMesh.pdf
endif

xi = linspace(0.2,1.1,5); yi = xi*0.8+0.05;
Ngrid = 100;
[xi,yi] = meshgrid(linspace(0,1,Ngrid),linspace(0,2,Ngrid));

function res = u_exact(xy)
   res = exp(xy(:,2));
endfunction

function u = f(xy)
  u = -exp(xy(:,2));
endfunction
u_ex = reshape(u_exact([xi(:),yi(:)]),Ngrid,Ngrid);

u = BVP2Dsym(Mesh,1,0,'f','u_exact',0,0);
[ui,uxi,uyi] = FEMgriddata(Mesh,u,xi,yi);

figure(2)
FEMtrimesh(Mesh,u)
hold on
plot3(xi,yi,ui,'g.')
hold off
xlabel('x'); ylabel('y'); title('u')
view([-60 25])

if printing
%%  print -dpdfwrite doc/InsideElementSolution.pdf
  print -dpng doc/InsideElementSolution.png
endif

figure(3)
mesh(xi,yi,ui-u_ex)
xlabel('x'); ylabel('y'); zlabel('difference u')
view([-100 25])

if printing
  switch MeshType
    case 'linear'
      %%    print -dpdfwrite doc/InsideElementDiffLin.pdf
    print -dpng doc/InsideElementDiffLin.png
  case 'quadratic'
    %%    print -dpdfwrite doc/InsideElementDiffQuad.pdf
    print -dpng doc/InsideElementDiffQuad.png
  case 'cubic'
    %%    print -dpdfwrite doc/InsideElementDiffCubic.pdf
    print -dpng doc/InsideElementDiffCubic.png
  endswitch
endif

figure(12)
mesh(xi,yi,uxi)
xlabel('x'); ylabel('y'); title('u_x')

figure(13)
mesh(xi,yi,uyi)
xlabel('x'); ylabel('y'); title('u_y')

if printing
  switch MeshType
    case 'linear'
      %%    print -dpdfwrite doc/InsideElementUyLin.pdf
      print -dpng doc/InsideElementUyLin.png
    case 'quadratic'
      %%    print -dpdfwrite doc/InsideElementUyQuad.pdf
      print -dpng doc/InsideElementUyQuad.png
    case 'cubic'
      %%    print -dpdfwrite doc/InsideElementUyCubic.pdf
      print -dpng doc/InsideElementUyCubic.png
  endswitch
endif

figure(22)
plot3(xi,yi,uxi,'k.')
xlabel('x'); ylabel('y'); title('u_x')

figure(23)
plot3(xi,yi,uyi,'k.')
xlabel('x'); ylabel('y'); title('u_y')

figure(14)
mesh(xi,yi,uyi-u_ex)
xlabel('x'); ylabel('y'); title('difference of u_y')
view([-110, 30])
if printing
  switch MeshType
    case 'linear'
      %%    print -dpdfwrite doc/InsideElementUyLin.pdf
      print -dpng doc/InsideElementUyDiffLin.png
    case 'quadratic'
      %%    print -dpdfwrite doc/InsideElementUyQuad.pdf
      print -dpng doc/InsideElementUyDiffQuad.png
    case 'cubic'
      %%    print -dpdfwrite doc/InsideElementUyCubic.pdf
      print -dpng doc/InsideElementUyDiffCubic.png
  endswitch
endif
