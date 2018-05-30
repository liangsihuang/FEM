function enf = BernoulliBeam2D_ElementNodeForce( ie )
% ����:
% ie:�ڵ��
% �����
% enf����Ԫ�ֲ�����ϵ�µĽڵ��� 
global Element dU DF
i = Element( ie, 1 ) ;
j = Element( ie, 2 ) ;
% ��ȡ��Ԫλ�ƣ�ע�⣺λ�����������꣡
d_element = zeros( 6, 1 ) ;
d_element( 1:3 ) = dU( (i-1)*3+1:(i-1)*3+3 ) ;
d_element( 4:6 ) = dU( (j-1)*3+1:(j-1)*3+3 ) ;
k = BernoulliBeam2D_Stiffness(ie,1);%1�����������꣬��Ϊλ�������������
enf = k * d_element ;
% ��Ҫ��ȥ��Ч�ڵ���
[num,~] = size(DF);
for idf = 1:1:num
    if ie == DF( idf, 1 ) 
        enf = enf - BernoulliBeam2D_EquivalentNodeForce( DF(idf,1), ...
              DF(idf, 2), DF( idf, 3), DF( idf, 4 ) ) ;
        break ;
    end
end
% �ѵ�Ԫ�ڵ���ת�����ֲ�������     
T = Beam2D_TransformMatrix( ie ) ;
enf = T  * enf ;
return