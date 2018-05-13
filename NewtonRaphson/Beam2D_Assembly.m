function gK=Beam2D_Assembly( ie,k,gK )
% �ѵ�Ԫ�նȾ��󼯳ɵ�����նȾ���
% ����:
% ie�� ��Ԫ��
% k����Ԫ�նȾ���
% �����
% gK������նȾ���
global gElement
for i=1:1:2
    for j=1:1:2
        for p=1:1:3
            for q =1:1:3
                m = (i-1)*3+p ;
                n = (j-1)*3+q ;
                M = (gElement(ie,i)-1)*3+p ;
                N = (gElement(ie,j)-1)*3+q ;
                gK(M,N) = gK(M,N) + k(m,n) ;
            end
        end
    end
end
return