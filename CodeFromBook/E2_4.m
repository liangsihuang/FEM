%
% Example 2.4 Divergence of the Newton-Raphson method
%
xdata=zeros(40,1);
ydata=zeros(40,1);
tol = 1.0e-5;
iter = 0;
u = .4990807536;  %This initial value will make oscillate
u = 0.5;
uold = u;
c=0;
P = u+atan(5*u);
R = -P;
conv= R^2;
xdata(1)=u;
ydata(1)=P;
while conv > tol && iter < 20
 Kt = 1+5*(cos(atan(5*u)))^2;
 delu = R/Kt;
 u = uold + delu;
 P = u+atan(5*u);
 R = -P;
 conv= R^2;
 uold = u;
 iter = iter + 1;
 xdata(2*iter)=u;   ydata(2*iter)=0;
 xdata(2*iter+1)=u; ydata(2*iter+1)=P;
end
%
%
plot(xdata,ydata);
hold on;
x=[-1:0.1:1];
y=x+atan(5*x);
plot(x,y)
