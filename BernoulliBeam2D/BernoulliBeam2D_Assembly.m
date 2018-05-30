function BernoulliBeam2D_Assembly( ie, k )
%  把单元刚度矩阵集成到整体刚度矩阵
%  输入参数:
%      ie  --- 单元号
%      k   --- 单元刚度矩阵
%  返回值:
%      无
global Element K
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