function K=beam2d_stiffness530(E,A,I,L,cs,Ele_F1)
% 输入：E,A,I,L,cs,Ele_F1
% 输出：
% K: global stiffness matrix
F=Ele_F1(1,4);
M1=Ele_F1(1,3);
M2=Ele_F1(1,6);
T=[cs(1,1),cs(1,2),0,0,0,0;
    -cs(1,2),cs(1,1),0,0,0,0;
    0,0,1,0,0,0,0;
    0,0,0,cs(1,1),cs(1,2),0;
    0,0,0,-cs(1,2),cs(1,1),0;
    0,0,0,0,0,0,1];
% 弹性刚度矩阵
KE = [  E*A/L           0          0 -E*A/L           0          0
               0  12*E*I/L^3  6*E*I/L^2      0 -12*E*I/L^3  6*E*I/L^2
               0   6*E*I/L^2    4*E*I/L      0  -6*E*I/L^2    2*E*I/L
          -E*A/L           0          0  E*A/L           0          0
               0 -12*E*I/L^3 -6*E*I/L^2      0  12*E*I/L^3 -6*E*I/L^2
               0   6*E*I/L^2    2*E*I/L      0  -6*E*I/L^2    4*E*I/L] ;
% 几何刚度矩阵
KG = [  F/L                        0                 -M1/L   -F/L                        0                -M2/L
          0             12*F*I/A/L^3      6*F*I/A/L^2+F/10      0   -(12*F*I/A/L^3+6*F/5/L)    6*F*I/A/L^2+F/10
      -M1/L         6*F*I/A/L^2+F/10    4*F*I/A/L+2*F*L/15   M1/L       -(6*F*I/A/L^2+F/10)    2*F*I/A/L-F*L/30
       -F/L                        0                  M1/L    F/L                        0                 M2/L
          0   -(12*F*I/A/L^3+6*F/5/L)   -(6*F*I/A/L^2+F/10)     0      12*F*I/A/L^3+6*F/5/L  -(6*F*I/A/L^2+F/10)
      -M2/L         6*F*I/A/L^2+F/10      2*F*I/A/L-F*L/30   M2/L        -(6*F*I/A/L^2+F/10)   4*F*I/A/L+2*F*L/15];
K=T'*(KE+KG)*T;
end