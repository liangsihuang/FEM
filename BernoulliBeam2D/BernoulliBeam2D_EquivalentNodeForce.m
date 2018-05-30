function enf = BernoulliBeam2D_EquivalentNodeForce( ie, p1, p2, force_type )
% ����:
% ie:��Ԫ��
% p1:��һ���ڵ��ϵķֲ�������ֵ
% p2:�ڶ����ڵ��ϵķֲ�������ֵ
% force_type:�ֲ���������:1����ֲ���������2����ֲ���������3����ֲ����
% ���:
% enf:��������ϵ�µ�Ч�ڵ�������
global Element Node
enf = zeros( 6, 1 ) ; % ���� 6x1 �ĵ�Ч�ڵ�������
xi = Node( Element( ie, 1 ), 1 ) ;
yi = Node( Element( ie, 1 ), 2 ) ;
xj = Node( Element( ie, 2 ), 1 ) ;
yj = Node( Element( ie, 2 ), 2 ) ;
L = sqrt( (xj-xi)^2 + (yj-yi)^2 ) ;
if force_type==1 
    % �ֲ������� 
    enf( 1 ) = (2*p1+p2)*L/6 ;
    enf( 4 ) = (p1+2*p2)*L/6 ;
end
if force_type==2
    % �ֲ�������
    enf( 2 ) = (7*p1+3*p2)*L/20 ;
    enf( 3 ) = (3*p1+2*p2)*L^2/60 ;
    enf( 5 ) = (3*p1+7*p2)*L/20 ;
    enf( 6 ) = -(2*p1+3*p2)*L^2/60 ;
end
if force_type==3
    % �ֲ����
    enf( 2 ) = -(p1+p2)/2 ;
    enf( 3 ) = (p1-p2)*L/12 ;
    enf( 5 ) = (p1+p2)/2 ;
    enf( 6 ) = -(p1-p2)*L/12 ;
end
%  �ѵ�Ч�ڵ���ת��������������
T = Beam2D_TransformMatrix( ie );
enf = T' * enf;
return