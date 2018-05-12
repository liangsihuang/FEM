%
% Example 3.8 Uniaxial bar--total Lagrangian formulation
%
tol = 1.0e-5;  iter = 0;
u = 0; uold = u; f = 100; E = 200;
strain = u + 0.5*u^2;
stress = E*strain;
P = stress*(1+u);
R = f - P;
conv= R^2/(1+f^2);
fprintf('\n iter      u1      E11      S11         conv');
fprintf('\n %3d  %7.5f  %7.5f %8.3f %12.3e %7.5f',iter,u,strain,stress,conv);
while conv > tol && iter < 20
 Kt = E*(1+u)^2 + stress;
 delu = R/Kt;
 u = uold + delu;
 strain = u + 0.5*u^2;
 stress = E*strain;
 P = stress*(1+u);
 R = f - P;
 conv= R^2/(1+f^2);
 uold = u;
 iter = iter + 1;
 fprintf('\n %3d  %7.5f  %7.5f %8.3f %12.3e %7.5f',iter,u,strain,stress,conv);
end
