%
% Example 2.7 Two nonlinear springs (Secant method)
%
tol = 1.0e-5;    iter = 0;   c = 0;
u = [0.1; 0.1];  uold = u;
f = [0; 100];
P = [300*u(1)^2+400*u(1)*u(2)-200*u(2)^2+150*u(1)-100*u(2)
     200*u(1)^2-400*u(1)*u(2)+200*u(2)^2-100*u(1)+100*u(2)];
R = P - f;    Rold = R;
conv= (R(1)^2+R(2)^2)/(1+f(1)^2+f(2)^2);
fprintf('\n iter      u1      u2         conv       c');
fprintf('\n %3d  %7.5f %7.5f %12.3e %7.5f',iter,u(1),u(2),conv,c);
Kt = [600*u(1)+400*u(2)+150  -400*u(2)+400*u(1)-100
      400*u(1)-400*u(2)-100  400*u(2)-400*u(1)+100];
while conv > tol && iter < 20
 delu = -Kt\R;
 u = uold + delu;
 P = [300*u(1)^2+400*u(1)*u(2)-200*u(2)^2+150*u(1)-100*u(2); 
      200*u(1)^2-400*u(1)*u(2)+200*u(2)^2-100*u(1)+100*u(2)];
 R = P - f;
 conv= (R(1)^2+R(2)^2)/(1+f(1)^2+f(2)^2);
 c = abs(0.9-u(2))/abs(0.9-uold(2))^2;
 delR = R - Rold;
 Kt = Kt + (delR-Kt*delu)*delu'/norm(delu)^2;
 uold = u;  Rold = R;
 iter = iter + 1;
 fprintf('\n %3d  %7.5f %7.5f %12.3e %7.5f',iter,u(1),u(2),conv,c);
end