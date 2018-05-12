%
% Example 2.6 Nonlinear algebraic equation (secant method)
%
xdata=zeros(40,1); ydata=zeros(40,1);
tol = 1.0e-5;    iter = 0;   c = 0;
u = 0.5;         uold = u;
P = u+atan(5*u); Pold = P;
R = -P;          conv= R^2;
fprintf('\n iter       u         conv       c');
fprintf('\n %3d  %7.5f %12.3e %7.5f',iter,u,conv,c);
Ks = 1+5*(cos(atan(5*u)))^2;
xdata(1)=u;
ydata(1)=P;
while conv > tol && iter < 20
 delu = R/Ks;
 u = uold + delu;
 P = u+atan(5*u);
 R = -P;
 conv= R^2;
 c = abs(u)/abs(uold)^2;
 Ks = (P - Pold)/(u - uold);
 uold = u;
 Pold = P;
 iter = iter + 1;
 xdata(2*iter)=u;   ydata(2*iter)=0;
 xdata(2*iter+1)=u; ydata(2*iter+1)=P;
 fprintf('\n %3d  %7.5f %12.3e %7.5f',iter,u,conv,c);
end
plot(xdata,ydata);
hold on;
x=[-2:0.1:2];
y=x+atan(5*x);
plot(x,y)
