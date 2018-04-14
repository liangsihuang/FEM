function enf = BernoulliBeam2D_ElementNodeForce( ie )
%  ���㵥Ԫ�Ľڵ���
%  �������
%      ie  ----- �ڵ��
%  ����ֵ
%      enf ----- ��Ԫ�ֲ�����ϵ�µĽڵ��� 
    global gElement gNode gDelta gDF
    i = gElement( ie, 1 ) ;
    j = gElement( ie, 2 ) ;
    de = zeros( 6, 1 ) ;
    de( 1:3 ) = gDelta( (i-1)*3+1:(i-1)*3+3 ) ;
    de( 4:6 ) = gDelta( (j-1)*3+1:(j-1)*3+3 ) ;
    k = BernoulliBeam2D_Stiffness(ie,1);%1������������
    enf = k * de ;
%     ��Ҫ��ȥ��Ч�ڵ���
    [df_number, dummy] = size( gDF ) ;
    for idf = 1:1:df_number
        if ie == gDF( idf, 1 ) 
            enf = enf - BernoulliBeam2D_EquivalentNodeForce( gDF(idf,1), ...
                  gDF(idf, 2), gDF( idf, 3), gDF( idf, 4 ) ) ;
            break ;
        end
    end
    
    T = BernoulliBeam2D_TransformMatrix( ie ) ;
    enf = transpose( T ) * enf ;
return