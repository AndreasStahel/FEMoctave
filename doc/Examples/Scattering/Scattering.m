## -*- texinfo -*-
## @deftypefn  {} {} Scattering.m
##
## This is a demo file inside the `doc/Examples/Scattering/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

%% a scattering example
global k len
k = 2*pi; d = 0.5; len = 3; a = d/10;
xy = [len d/2 -2; len -d/2 -2; -len -d/2 -2; -len d/2 -2];
Np = 36;
phi = linspace(0,2*pi,Np+1)(1:Np); phi = phi';
Wire.name = 'Hole';
Wire.border = [a*cos(phi) a*sin(phi) -ones(Np,1)];
Wire.point = [0 0];
Mesh = CreateMeshTriangle('strip',xy,1/(200),Wire);
Mesh = MeshUpgrade(Mesh,'cubic');

function res = gn1(xy,dummy)
   global k len
   res = 2*j*k*(xy(:,1)<=(-len+2*eps));      %% applies on left edge
end
function res = gn2(xy,dummy)
   global k len
   res = -j*k*(abs(xy(:,1))>=(len-2*eps));  %% applies on left and right edge
end

u = BVP2D(Mesh,1,-k^2,0,0,0,0,'gn1','gn2','type','complex');

figure(1); FEMtrisurf(Mesh,real(u)); ylim([-d/2,d/2])
           xlabel('x'); ylabel('y'); zlabel('real(u(x))')
figure(2); FEMtrisurf(Mesh,imag(u)); ylim([-d/2,d/2])
           xlabel('x'); ylabel('y'); zlabel('imag(u(x))')
figure(3); FEMtrisurf(Mesh,abs(u)); ylim([-d/2,d/2])
           xlabel('x'); ylabel('y'); zlabel('abs(u(x))')

x = linspace(-len,+len,5001); y = zeros(size(x));
[u_mid] = FEMgriddata(Mesh,u,x,y);
figure(11); plot(x,real(u_mid)); xlabel('x');  title('real(u(x))')
figure(12); plot(x,imag(u_mid)); title('imag(u(x))')
figure(13); plot(x,abs(u_mid)); title('abs(u(x))')
figure(14); plot3(x,real(u_mid),imag(u_mid));
            xlabel('x'); ylabel('real(u)'); zlabel('imag(u)'); axis equal
figure(15); plot(real(u_mid),imag(u_mid));
            xlabel('real(u)'); ylabel('imag(u)'); axis equal

tau_transmit = u_mid(end)
tau_reflect  = u_mid(1)-1
