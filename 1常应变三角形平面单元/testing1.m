clear;clc;
% % 节点：节点号 x坐标 y坐标
% Node=[1 2.0 1.0;
%     2 2.0 0.0;
%     3 0.0 1.0;
%     4 0.0 0.0];
% % 单元：单元号 3个节点号（逆时针） ID 材料号
% element=[1 2 3 4;
%     2 2 1 3];
% 单位：N，m
E=1e7;
NU=1/3;
t=0.1;
ID=1;
k1=Triangle2D3Node_Stiffness(E,NU,t,2,0,0,1,0,0,ID);
k2=Triangle2D3Node_Stiffness(E,NU,t,0,1,2,0,2,1,ID);
KK=zeros(8,8);
KK=Triangle2D3Node_Assembly(KK,k1,2,3,4);
KK=Triangle2D3Node_Assembly(KK,k2,3,2,1);
P=zeros(8,1);
P(2,1)=-50000;
P(4,1)=-50000;
[KK,P]=Triangle2D3Node_Boundary(KK,P,5,0);
[KK,P]=Triangle2D3Node_Boundary(KK,P,6,0);
[KK,P]=Triangle2D3Node_Boundary(KK,P,7,0);
[KK,P]=Triangle2D3Node_Boundary(KK,P,8,0);
u=KK\P;
