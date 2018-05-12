function enf = NonlinearBeam2D_ElementNodeForce( ie )
% ����:
% ie:�ڵ��
% �����
% enf����Ԫ�ֲ�����ϵ�µĽڵ��� 
global gElement gDelta
i = gElement( ie, 1 ) ;
j = gElement( ie, 2 ) ;
% ��ȡ��Ԫλ�ƣ�ע�⣺λ�����������꣡
d_element = zeros( 6, 1 ) ;
d_element( 1:3 ) = gDelta( (i-1)*3+1:(i-1)*3+3 ) ;
d_element( 4:6 ) = gDelta( (j-1)*3+1:(j-1)*3+3 ) ;
k = NonlinearBeam2D_Stiffness(ie,1);%1�����������꣬��Ϊλ�������������
enf = k * d_element;
% �ȼ���û�е�Ч�ڵ���......

% �ѵ�Ԫ�ڵ���ת�����ֲ�������     
T = Beam2D_TransformMatrix( ie ) ;
enf = transpose( T ) * enf ;
return