function F=InternalForce(du)
% ���룺
% du: ���������µ�λ���������󣬽ڵ���*3
% �����
% F: ���������µĽڵ��������������󣬽ڵ���*3
global Node Element InF

% ��ȡ��Ԫ��du
[num,~]=sizes(Element);
for i=1:num
    node_elem(1:2)=Node(n1*2-1:n1*2);
    node_elem(3:4)=Node(n2*2-1:n2*2);
    du_purified=RemoveRigidMotion(node_elem,du_elem);
    node_elem_C1=node_elem-du_elem(1:2,4:5);
    inf=zeros(6,1);
    k=NonlinearBeam2D_Stiffness(ie,2,inf);
    
end