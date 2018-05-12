%
% Example 4.15 - Shear deformation of elastoplastic square
%
Young = 24000; nu=0.2; mu=Young/2/(1+nu); lambda=nu*Young/((1+nu)*(1-2*nu));
beta = 0; H = 1000; sY = 200*sqrt(3);
mp = [lambda mu beta H sY];
Iden=[1 1 1 0 0 0]';
coef=Young/((1+nu)*(1-2*nu));
D=2*mu*eye(6) + lambda*Iden*Iden';
D(4,4) = mu; D(5,5) = mu; D(6,6) = mu;
stressN=[0 0 0 0 0 0]';
deps=[0 0 0 0 0 0]';
alphaN = [0 0 0 0 0 0]';
epN=0;
for i=1:15
    deps(4) = 0.004;
    [stress, alpha, ep]=combHard(mp,D,deps,stressN,alphaN,epN);
    X(i) = i*deps(4); Y(i) = stress(4); Z(i) = ep;
    stressN = stress; alphaN = alpha; epN = ep;
end
X = [0 X]; Y=[0 Y]; plot(X,Y);
