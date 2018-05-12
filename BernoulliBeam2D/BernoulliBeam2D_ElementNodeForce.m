function enf = BernoulliBeam2D_ElementNodeForce( ie )
% ����:
% ie:�ڵ��
% �����
% enf����Ԫ�ֲ�����ϵ�µĽڵ��� 
global gElement gDelta gDF
i = gElement( ie, 1 ) ;
j = gElement( ie, 2 ) ;
% ��ȡ��Ԫλ�ƣ�ע�⣺λ�����������꣡
d_element = zeros( 6, 1 ) ;
d_element( 1:3 ) = gDelta( (i-1)*3+1:(i-1)*3+3 ) ;
d_element( 4:6 ) = gDelta( (j-1)*3+1:(j-1)*3+3 ) ;
k = BernoulliBeam2D_Stiffness(ie,1);%1�����������꣬��Ϊλ�������������
enf = k * d_element ;
% ��Ҫ��ȥ��Ч�ڵ���
[distributedforce_number,~] = size(gDF);
for idf = 1:1:distributedforce_number
    if ie == gDF( idf, 1 ) 
        enf = enf - BernoulliBeam2D_EquivalentNodeForce( gDF(idf,1), ...
              gDF(idf, 2), gDF( idf, 3), gDF( idf, 4 ) ) ;
        break ;
    end
end
% �ѵ�Ԫ�ڵ���ת�����ֲ�������     
T = BernoulliBeam2D_TransformMatrix( ie ) ;
enf = transpose( T ) * enf ;
return