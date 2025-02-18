clear *
global k len a
k = 2*pi; d = 0.5; len = 4; a = d/3;
xy = [len d/2 -2; len -d/2 -2; -len -d/2 -2; -len d/2 -2];
Mesh = CreateMeshTriangle('strip',xy,len/(2*800));
Mesh = MeshUpgrade(Mesh,'cubic');

function res = gn1(xy,dummy)
   global k len
   res = 2*j*k*(xy(:,1)<=(-len+10*eps));         %% applies on left edge
end

function res = gn2(xy,dummy)
   global k len
   res = (0-1j*k)*(abs(xy(:,1))>=(len-10*eps));  %% applies on left and right edge
end

function res = k2(xy)
  global k a
  alpha = 0.25;
  x = xy(:,1); y = xy(:,2);
  res = -k^2*ones(size(x));
  if 1  %% circle with different speed
     ind = (x.^2+y.^2)<a^2;
   else  %% band with different speed
     ind = x.^2<a^2;
  endif
  res(ind) *= alpha;
endfunction

u = BVP2D(Mesh,1,'k2',0,0,0,0,'gn1','gn2','type','complex');

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

