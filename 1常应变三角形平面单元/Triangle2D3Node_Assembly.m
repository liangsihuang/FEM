function z=Triangle2D3Node_Assembly(KK,k,i,j,m)
% 该函数进行单元刚度矩阵的组装
% 输入单元刚度矩阵k
% 输入单元的节点编号i,j,m
% 输出整体刚度矩阵KK

% 计算送入整体刚度矩阵中的位置，2为每个节点两个自由度
DOF(1)=2*i-1;
DOF(2)=2*i;
DOF(3)=2*j-1;
DOF(4)=2*j;
DOF(5)=2*m-1;
DOF(6)=2*m;
% 送入整体刚度矩阵
for n1=1:6
    for n2=1:6
        KK(DOF(n1),DOF(n2))=KK(DOF(n1),DOF(n2))+k(n1,n2);
    end
end
z=KK;
end