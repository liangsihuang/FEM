function enf = NonlinearBeam2D_ElementNodeForce( ie )
% ����:
% ie:�ڵ��
% �����
% enf����Ԫ�ֲ�����ϵ�µĽڵ���
deform=RemoveRigidMotion(ie);
inf=zeros(6,1); %�����Ӧ��Ϊ0
k = NonlinearBeam2D_Stiffness(ie,2,inf);%2����ֲ�����
enf = k * deform;
return