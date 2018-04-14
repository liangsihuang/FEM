function k=BernoulliBeam2D_Stiffness(ie,icoord)
%  ���㵥Ԫ�նȾ���
%  �������:
%     ie -------  ��Ԫ��
%     icoord  --  ����ϵ��������������������֮һ
%                    1  ----  ��������ϵ
%                    2  ----  �ֲ�����ϵ
%  ����ֵ:
%     k  ----  ����icoord��ֵ����Ӧ����ϵ�µĸնȾ���
    global gNode gElement gMaterial
    k = zeros( 6, 6 ) ;
    E = gMaterial( gElement(ie, 3), 1 ) ;
    I = gMaterial( gElement(ie, 3), 2 ) ;
    A = gMaterial( gElement(ie, 3), 3 ) ;
    xi = gNode( gElement( ie, 1 ), 1 ) ;
    yi = gNode( gElement( ie, 1 ), 2 ) ;
    xj = gNode( gElement( ie, 2 ), 1 ) ;
    yj = gNode( gElement( ie, 2 ), 2 ) ;
    L = ( (xj-xi)^2 + (yj-yi)^2 )^(1/2) ;
    k = [  E*A/L           0          0 -E*A/L           0          0
               0  12*E*I/L^3  6*E*I/L^2      0 -12*E*I/L^3  6*E*I/L^2
               0   6*E*I/L^2    4*E*I/L      0  -6*E*I/L^2    2*E*I/L
          -E*A/L           0          0  E*A/L           0          0
               0 -12*E*I/L^3 -6*E*I/L^2      0  12*E*I/L^3 -6*E*I/L^2
               0   6*E*I/L^2    2*E*I/L      0  -6*E*I/L^2    4*E*I/L] ;
    if icoord == 1
        T = BernoulliBeam2D_TransformMatrix( ie ) ;
        k = T*k*transpose(T) ;
    end
return

end