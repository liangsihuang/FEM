function Triangle2D3Node_SolveModel
% step1. 定义整体刚度矩阵和节点力向量
[node_number,dummy] = size( gNode ) ;
KK=zeros(node_number*2,node_number*2);
P=zeros(node_number*2,1);
% step2. 计算单元刚度矩阵，并集成到整体刚度矩阵中
[element_number,dummy] = size( gElement ) ;
for ie=1:1:element_number
    k = StiffnessMatrix( ie, 1 ) ;
    AssembleStiffnessMatrix( ie, k ) ;
end
end