function T = BernoulliBeam2D_TransformMatrix( ie )
%  ���㵥Ԫ������ת������( �ֲ����� -> �������� )
%  �������
%      ie  ----- �ڵ��
%  ����ֵ
%      T ------- �Ӿֲ����굽�������������ת������
    global gElement gNode
    xi = gNode( gElement( ie, 1 ), 1 ) ;
    yi = gNode( gElement( ie, 1 ), 2 ) ;
    xj = gNode( gElement( ie, 2 ), 1 ) ;
    yj = gNode( gElement( ie, 2 ), 2 ) ;
    L = sqrt( (xj-xi)^2 + (yj-yi)^2 ) ;
    c = (xj-xi)/L ;
    s = (yj-yi)/L ;
    T=[ c  -s   0   0   0   0
        s   c   0   0   0   0
        0   0   1   0   0   0
        0   0   0   c  -s   0
        0   0   0   s   c   0
        0   0   0   0   0   1] ;
return