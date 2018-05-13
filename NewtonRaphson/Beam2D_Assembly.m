function gK=Beam2D_Assembly( ie,k,gK )
% 把单元刚度矩阵集成到整体刚度矩阵
% 输入:
% ie： 单元号
% k：单元刚度矩阵
% 输出：
% gK：整体刚度矩阵
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