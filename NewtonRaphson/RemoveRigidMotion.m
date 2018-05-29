function deform=RemoveRigidMotion(ie)
% ���룺
% ie����Ԫ��
% �����
% deform: C2�ֲ������µĵ��Ա��� (6,1)
global Element dU Node
n1=Element(ie,1);
n2=Element(ie,2);
du_elem(1:3)=dU(n1*3-2:n1*3);
du_elem(4:6)=dU(n2*3-2:n2*3);
T=Beam2D_TransformMatrix(ie);
du_elem=T'*du_elem; %ת������Ԫ�ֲ����꣬C1״̬
% C1��C2��������������
ua=du_elem(1);
wa=du_elem(2);
ra=du_elem(3);
ub=du_elem(4);
wb=du_elem(5);
rb=du_elem(6);
% ��C1״̬�¸˼��ĳ���, Node�����ǻ�δ����dU�Ľڵ�����
x1=Node(n1,1);
y1=Node(n1,2);
x2=Node(n2,1);
y2=Node(n2,2);
L1=sqrt((x1-x2)^2+(y1-y2)^2);
% C2״̬�³��ȣ�C2���C1״̬�ĸ���ת��
dX=L1+ub-ua;
dZ=wb-wa;
L2=sqrt(dX^2+dZ^2);
r=atan(dZ/dX);
% C2�ֲ������µĵ��Ա���
deform=zeros(6,1);
deform(4)=L2-L1;
deform(3)=ra-r;
deform(6)=rb-r;
end

