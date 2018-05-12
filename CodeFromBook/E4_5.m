%
% Example 4.5 Two bars in parallel: updating iteration
%
E1=10000; Et1=1000; sYield1=5;
E2=5000;  Et2=500;  sYield2=7.5;
mp1 = [E1, 1, E1*Et1/(E1-Et1), sYield1];
mp2 = [E2, 0, E2*Et2/(E2-Et2), sYield2];
nS1 = 0; nA1 = 0; nep1 = 0;
nS2 = 0; nA2 = 0; nep2 = 0;
A1 = 0.75; L1 = 100;
A2 = 1.25; L2 = 100;
tol = 1.0E-5; u = 0; P = 15; iter = 0;
Res = P - nS1*A1 - nS2*A2;
Dep1 = E1; Dep2 = E2;
conv = Res^2/(1+P^2);
fprintf('\niter        u      S1      S2      A1      A2');
fprintf('      ep1      ep2   Residual');
fprintf('\n %3d  %7.4f %7.3f %7.3f %7.3f %7.3f %8.6f %8.6f %10.3e',...
      iter,u,nS1,nS2,nA1,nA2,nep1,nep2,Res);
while conv > tol && iter < 20
  delu = Res / (Dep1*A1/L1 + Dep2*A2/L2);
  u = u + delu;
  delE = delu / L1;
  [Snew1, Anew1, epnew1]=combHard1D(mp1,delE,nS1,nA1,nep1);
  [Snew2, Anew2, epnew2]=combHard1D(mp2,delE,nS2,nA2,nep2);
  Res = P - Snew1*A1 - Snew2*A2;
  conv = Res^2/(1+P^2);
  iter = iter + 1;
  Dep1 = E1; if epnew1 > nep1; Dep1 = Et1; end
  Dep2 = E2; if epnew2 > nep2; Dep2 = Et2; end
  nS1 = Snew1; nA1 = Anew1; nep1 = epnew1;
  nS2 = Snew2; nA2 = Anew2; nep2 = epnew2;
%  fprintf('\n %3d  %7.4f %7.3f %7.3f %7.3f %7.3f %8.6f %8.6f %10.3e',...
  fprintf('\n %3d  %8.5f %8.5f %8.5f %8.5f %8.5f %9.7f %9.7f %10.3e',...
      iter,u,nS1,nS2,nA1,nA2,nep1,nep2,Res);
end
