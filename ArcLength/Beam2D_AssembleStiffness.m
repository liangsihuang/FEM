function K=Beam2D_AssembleStiffness(K,ie,k)
% 把单元刚度矩阵集成到整体刚度矩阵
% 输入:
% ie： 单元号
% k：单元刚度矩阵
% 输出：
% K：整体刚度矩阵
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