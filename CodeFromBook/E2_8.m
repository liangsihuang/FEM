%
% Example 2.8 Displacement controlled procedure
%
tol = 1.0e-5;  conv = 0;  u1 = 0;  u1old = u1;
fprintf('\n step      u1      u2       F');
% Displacement increment loop
for i=1:9
 u2 = 0.1*i;
 P = 300*u1^2+400*u1*u2-200*u2^2+150*u1-100*u2;
 R = -P;
 conv = R^2;
 % Convergence loop
 iter = 0;
 while conv > tol && iter < 20
  Kt = 600*u1+400*u2+150;
  delu1 = R/Kt;
  u1 = u1old + delu1;
  P = 300*u1^2+400*u1*u2-200*u2^2+150*u1-100*u2;
  R = -P;
  conv= R^2;
  u1old = u1;
  iter = iter + 1;
 end
 F = 200*u1^2-400*u1*u2+200*u2^2-100*u1+100*u2;
 fprintf('\n %3d  %7.5f %7.5f %7.3f',i,u1,u2,F);
end
