%
% Example 3.9 Uniaxial bar--updated Lagrangian formulation
%
tol = 1.0e-5;  iter = 0;  E = 200;
u = 0;         uold = u;  f = 100;
strain = u/(1+u);
stress = E*(u+.5*u^2)*(1+u);
P = stress;
R = f - P;
conv= R^2/(1+f^2);
fprintf('\n iter      u1      E11      S11         conv');
fprintf('\n %3d  %7.5f  %7.5f %8.3f %12.3e %7.5f',iter,u,strain,stress,conv);
while conv > tol && iter < 20
 Kt = E*(1+u)^2 + stress/(1+u);
 delu = R/Kt;
 u = uold + delu;
 strain = u/(1+u);
 stress = E*(u+.5*u^2)*(1+u);
 P = stress;
 R = f - P;
 conv= R^2/(1+f^2);
 uold = u;
 iter = iter + 1;
 fprintf('\n %3d  %7.5f  %7.5f %8.3f %12.3e %7.5f',iter,u,strain,stress,conv);
end
