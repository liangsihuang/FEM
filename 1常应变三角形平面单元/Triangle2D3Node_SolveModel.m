function Triangle2D3Node_SolveModel
% step1. ��������նȾ���ͽڵ�������
[node_number,dummy] = size( gNode ) ;
KK=zeros(node_number*2,node_number*2);
P=zeros(node_number*2,1);
% step2. ���㵥Ԫ�նȾ��󣬲����ɵ�����նȾ�����
[element_number,dummy] = size( gElement ) ;
for ie=1:1:element_number
    k = StiffnessMatrix( ie, 1 ) ;
    AssembleStiffnessMatrix( ie, k ) ;
end
end