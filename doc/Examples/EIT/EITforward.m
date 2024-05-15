## -*- texinfo -*-
## @deftypefn  {} {} EITforward.m
##
## This is a demo file  inside the `doc/Examples/EIT/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

global Rx Ry dalpha my_angle
N = 2*64; %% number of angle segments
alpha = linspace(0,2*pi*(N-1)/N,N)';  Rx = 1; Ry = 0.5;
dalpha = 2*(alpha(2)-alpha(1));
x = Rx*cos(alpha); y = Ry*sin(alpha);
BC = -2*ones(size(x));
my_angle = 120  %% select the configuration, use 60 or 120

function res = sigma(xy)                 %% the conductivity
  x = xy(:,1); y = xy(:,2);
  res = ones(size(x));
  res((x+0.5).^2+y.^2<=0.25^2) *=  4 ;   %% heart on the left
  res((x-0.4).^2+y.^2<=0.35^2) *=  1/4;  %% lung on the right
endfunction

FEMmesh = CreateMeshTriangle('EIT',[x,y,BC],0.003);
FEMmesh = MeshUpgrade(FEMmesh,'cubic');

figure(1); FEMtrimesh(FEMmesh,sigma(FEMmesh.nodes));  %% show the conductivity
           xlabel('x'); ylabel('y'); zlabel('\sigma'); view(20,50)

function res = flux_n(xy)         %% define the current density on the boundary
  global dalpha my_angle Rx Ry
  alpha = atan2(xy(:,2)/Ry,xy(:,1)/Rx);  %% assure correct angle
  res = zeros(size(alpha));
  res(abs(alpha+pi/2) < dalpha) = -1;
  switch my_angle
    case 60
      res(abs(alpha-pi/3)   < dalpha) = +1;
    case 120
      res(abs(alpha-pi*2/3) < dalpha) = +1;
  endswitch
  res = res./sqrt(Rx^2*sin(alpha).^2 + Ry^2*cos(alpha).^2);  %% adjust for the arc length
endfunction

u_0 = BVP2DsymMean(FEMmesh,   1   ,0,0,0,'flux_n',0);  %% the reference result
u   = BVP2DsymMean(FEMmesh,'sigma',0,0,0,'flux_n',0);  %% the actual result

figure(2); FEMtrimesh(FEMmesh,u)                    %% show the solution
           xlabel('x'); ylabel('y');

figure(3); clf; FEMtricontour(FEMmesh,u,41)         %% show the contour levels
           hold on;
           plot([x;x(1)],[y;y(1)],'k');             %% add the boundary
           hold off
           xlabel('x'); ylabel('y'); axis equal

u_boundary   = FEMgriddata(FEMmesh,u,  x,y);
u_0_boundary = FEMgriddata(FEMmesh,u_0,x,y);

figure(4); plot(alpha*180/pi,u_boundary,alpha*180/pi,u_0_boundary)
           xlabel('angle [deg]'); ylabel('u'); xlim([0,360])
           legend('true','reference')     %% show the voltages on the boundary

figure(5); plot(alpha*180/pi,u_boundary-u_0_boundary)
           xlabel('angle [deg]'); ylabel('u-u_0'); xlim([0,360])
           legend('true-reference')  %% show the difference

%% create the vector field for the current density
[xx,yy] = meshgrid(linspace(-Rx,Rx,21),linspace(-0.8*Ry,0.8*Ry,21));
[ui,uxi,uyi] = FEMgriddata(FEMmesh,u,xx,yy);
conductivity = reshape(sigma([xx(:),yy(:)]),size(xx));
uxi = conductivity.*uxi;
uyi = conductivity.*uyi;

figure(6); quiver(xx,yy,uxi,uyi,2)    %%% show the vector field
           xlabel('x'); ylabel('y')
           hold on;
           plot([x;x(1)],[y;y(1)],'k');  %% add the boundary
           hold off
%% create and show the streamlines
streamline(xx,yy,uxi,uyi,[-0.3 -0.2,-0.1,0,0.1,0.2 0.3],-0.8*Ry*ones(1,7));
