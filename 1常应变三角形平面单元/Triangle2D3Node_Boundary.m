function [z,p]=Triangle2D3Node_Boundary(KK,P,DOF,disp)
% 该函数用乘大数法修改总刚度矩阵，实现边界条件
% 输入总刚度矩阵KK，荷载向量P，要约束的自由度DOF，约束位移值disp
% 输出总刚度矩阵z，荷载向量p

BigNumber=10^20;
P(DOF,1)=BigNumber*KK(DOF,DOF)*disp;
p=P;
KK(DOF,DOF)=KK(DOF,DOF)*BigNumber;
z=KK;
end