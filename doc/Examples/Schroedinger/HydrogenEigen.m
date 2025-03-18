## -*- texinfo -*-
## @deftypefn  {} {} HydrogenEigen.m
##
## This is a demo file inside the `doc/Examples/Schroedinger/` directory@*
## Find the description in the documentation FEMdoc.pdf
##
## @end deftypefn

%% find eigenvalues of the hydrogen atom
n = 6; %% number of eigenvalues

L = 200;
pkg load femoctave
interval = linspace(0,L,1000);
shift = 2;
a = 1; b = 0; c = @(x)-1./x+shift; w = 1;
[x,EVal,EVec] = BVP1Deig(interval,a,b,c,w,0,0,n);

EVal = EVal-shift;
OneOverEvalues = 1./EVal

[MaxVal,MaxInd] = max(abs(EVec));
Sign = zeros(1,n);
for ii=1:n
  Sign(ii) = sign(EVec(MaxInd(ii),ii));
endfor

r = x/2; v1 = abs(EVec(:,1)./x);
figure(1); plot(r,EVec*diag(Sign))
           xlabel('r/a_0'); ylabel('v(r)'); xlim([0,30])
figure(2); plot(r,log(v1))
           xlabel('r/a_0'); ylabel('ln(v_1(r)/r)'); xlim([0,30])

printing = 0;
if printing
  figure(1); print -dpdfcrop HydrogenEigenGraph.pdf
  figure(2); print -dpdfcrop HydrogenEigenLog_u1.pdf
endif

