function K=Beam2D_AssembleStiffness(K,ie,k)
% �ѵ�Ԫ�նȾ��󼯳ɵ�����նȾ���
% ����:
% ie�� ��Ԫ��
% k����Ԫ�նȾ���
% �����
% K������նȾ���
global Element
for i=1:1:2
    for j=1:1:2
        for p=1:1:3
            for q =1:1:3
                m = (i-1)*3+p ;
                n = (j-1)*3+q ;
                M = (Element(ie,i)-1)*3+p ;
                N = (Element(ie,j)-1)*3+q ;
                K(M,N) = K(M,N) + k(m,n) ;
            end
        end
    end
end
return