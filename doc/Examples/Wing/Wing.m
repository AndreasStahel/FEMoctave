## -*- texinfo -*-
## @deftypefn  {} {} Wing.m
##
## This is a demo file inside the `doc/Examples/Wing/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

clear *
x = linspace(0,1,200)';
function [height,slope] = FoilThickness(x)
  FoilT = 0.10; %% thickness of profile
  height =  5*FoilT*(0.2969*sqrt(x)-0.1260*x-0.3516*x.^2+0.2843*x.^3- 0.1015*x.^4);
  fix = height(end)/(5*FoilT);
  height =  5*FoilT*(0.2969*sqrt(x)-0.1260*x-0.3516*x.^2+0.2843*x.^3- (0.1015+fix)*x.^4);
  slope =   5*FoilT*(0.2969*0.5./sqrt(x)-0.1260-2*0.3516*x+3*0.2843*x.^2- 4*(0.1015+fix)*x.^3);
  slope(1) = (2*sqrt(2)-1)*slope(2);
endfunction
[height,slope] = FoilThickness(x);
figure(1); plot(x,[-height,height]); axis equal
Angle = -10/180*pi;  %% angle of attack
R = [cos(Angle),-sin(Angle);+sin(Angle),cos(Angle)];
Upper = (R*[x,height]')';  Lower = (R*[x,-height]')';
AngleNormalUpper = atan2(+1,-slope)+Angle;
AngleNormalLower = atan2(-1,-slope)+Angle;

DomainHole = [flipud(Upper);Lower(1:end-1,:)];
Hole.name = 'hole';
Hole.border = [DomainHole,-2*ones(length(DomainHole),1)];
Hole.point = [0.1,0];
Borders = [-0.5,-0.6,-2;2,-0.6,-1;2,0.5,-2;-0.5,0.5,-1];
Mesh = CreateMeshTriangle('Wing',Borders,1e-2,Hole);
figure(2); FEMtrimesh(Mesh); axis equal
Mesh = MeshUpgrade(Mesh,'cubic');

function res = gD(xy)
  res = -xy(:,1);
endfunction

u = BVP2Dsym(Mesh,1,0,0,'gD',0,0);
[ux,uy] = FEMEvaluateGradient(Mesh,u);
v = sqrt(sum([ux.^2,uy.^2],2));
figure(3); FEMtrimesh(Mesh,v.^2); xlabel('x'); ylabel('y'); zlabel('v^2');
           zlim([0,2]); caxis([min(v.^2),2])
Levels = [0.5:0.1:2];
figure(4); clf; FEMtricontour(Mesh,v.^2,Levels); xlabel('x'); ylabel('y');
title('contours of v^2'); colorbar(); axis equal
           hold on; plot(DomainHole(:,1),DomainHole(:,2),'k'); hold off

x = x(1:end-1); Upper = Upper(1:end-1,:); Lower = Lower(1:end-1,:);
AngleNormalUpper = AngleNormalUpper(1:end-1);
AngleNormalLower = AngleNormalLower(1:end-1);

function res = ArcLength(x,y);
  dx = diff(x); dy = diff(y); ds = sqrt(dx.^2+dy.^2);
  res = [0;cumsum(ds)];
endfunction

dsUpper = ArcLength(Upper(:,1),Upper(:,2));
pUpper = FEMgriddata(Mesh,v.^2,Upper(:,1),Upper(:,2)).*sin(AngleNormalUpper);
ForceUpper = trapz(dsUpper,pUpper)
dsLower = ArcLength(Lower(:,1),Lower(:,2));
pLower = FEMgriddata(Mesh,v.^2,Lower(:,1),Lower(:,2)).*sin(AngleNormalLower);
ForceLower = trapz(dsLower,pLower)

figure(5); plot(dsUpper,pUpper,dsLower,-pLower); ylabel('pressure');
           legend('p_{upper}','-p_{lower}','location','north');
           xlim([0,1]); ylim([0,max(pUpper)])


[xx,yy] = meshgrid(linspace(-0.5,1.5,201),linspace(-0.5,0.5,201));
[ui,uxi,uyi] = FEMgriddata(Mesh,u,xx,yy);

figure(6); clf; NN = 25;
           streamline(xx,yy,-uxi,-uyi,-0.5*ones(1,NN),linspace(-0.3,0.2,NN),[0.1,10000]);
           hold on; plot(DomainHole(:,1),DomainHole(:,2),'k'); hold off; axis equal
