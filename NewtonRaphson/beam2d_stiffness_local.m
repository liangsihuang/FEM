% 子程序3：局部坐标刚度矩阵
function K=beam2d_stiffness_local(E,A,I,L,ELG)
F=ELG(1,4);
M1=ELG(1,3);
M2=ELG(1,6);
KL=[E*A/L,0,0,-E*A/L,0,0;... 
    0,12*E*I/L^3,6*E*I/L^2,0,-12*E*I/L^3,6*E*I/L^2;...
    0,6*E*I/L^2,4*E*I/L,0,-6*E*I/L^2,2*E*I/L;...
    -E*A/L,0,0,E*A/L,0,0;...
    0,-12*E*I/L^3,-6*E*I/L^2,0,12*E*I/L^3,-6*E*I/L^2;...
    0,6*E*I/L^2,2*E*I/L,0,-6*E*I/L^2,4*E*I/L];
KN=[F/L,0,-M1/L,-F/L,0,-M2/L;...
    0,12*F*I/A/L^3+6*F/5/L,6*F*I/A/L^2+F/10,0,-12*F*I/A/L^3-6*F/5/L,...
    6*F*I/A/L^2+F/10;...
    -M1/L,6*F*I/A/L^2+F/10,4*F*I/A/L+2*F*L/15,M1/L,-6*F*I/A/L^2-F/10,...
    2*F*I/A/L-F*L/30;...
    -F/L,0,M1/L,F/L,0,M2/L;...
    0,-12*F*I/A/L^3-6*F/5/L,-6*F*I/A/L^2-F/10,0,12*F*I/A/L^3+6*F/5/L,...
   -6*F*I/A/L^2-F/10;...
   -M2/L,6*F*I/A/L^2+F/10,2*F*I/A/L-F*L/30,M2/L,-6*F*I/A/L^2-F/10,...
   4*F*I/A/L+2*F*L/15];
K=KL+KN;
end