function Triangle2D3Node_SolveModel
global gNode gElement gMaterial gBC1 gNF gDF gK gDelta ID gVonMises
% step1. ��������նȾ���ͽڵ�������
[node_number,dummy] = size( gNode ) ;
gK=zeros(node_number*2,node_number*2);
P=zeros(node_number*2,1);
% step2. ���㵥Ԫ�նȾ��󣬲����ɵ�����նȾ�����
[element_number,dummy] = size( gElement ) ;
for ie=1:1:element_number
    k = Triangle2D3Node_Stiffness(ie,ID);
    gK= Triangle2D3Node_Assembly(gK,ie,k) ;
end
% step3. �Ѽ�����ֱ�Ӽ��ɵ�����ڵ���������
[nf_number, dummy] = size( gNF ) ;
for inf=1:1:nf_number
    n = gNF( inf, 1 ) ;
    d = gNF( inf, 2 ) ;
    P( (n-1)*2 + d ) = gNF( inf, 3 ) ;
end
% step5. ����Լ���������޸ĸնȾ���ͽڵ������������ó˴�����
[bc_number,dummy] = size( gBC1 ) ;
for ibc=1:1:bc_number
    n = gBC1(ibc, 1 ) ;
    d = gBC1(ibc, 2 ) ;
    m = (n-1)*2 + d ;
    P(m) = gBC1(ibc, 3)* gK(m,m) * 1e15 ;
    gK(m,m) = gK(m,m) * 1e15 ;
end
% step 6. ��ⷽ���飬�õ��ڵ�λ������
gDelta = gK \ P ;
% step 7. von-MisesӦ������
[element_number,dummy] = size( gElement ) ;
gVonMises=zeros(element_number,1);
for ie=1:1:element_number
    stress=Triangle2D3Node_Stress(ie,ID);
    gVonMises(ie,1)=Triangle2D3Node_vonMises(stress);
end
end