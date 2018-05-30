function enf = BernoulliBeam2D_ElementNodeForce( ie )
% 输入:
% ie:节点号
% 输出：
% enf：单元局部坐标系下的节点力 
global Element dU DF
i = Element( ie, 1 ) ;
j = Element( ie, 2 ) ;
% 提取单元位移，注意：位移是整体坐标！
d_element = zeros( 6, 1 ) ;
d_element( 1:3 ) = dU( (i-1)*3+1:(i-1)*3+3 ) ;
d_element( 4:6 ) = dU( (j-1)*3+1:(j-1)*3+3 ) ;
k = BernoulliBeam2D_Stiffness(ie,1);%1代表整体坐标，因为位移是整体坐标的
enf = k * d_element ;
% 还要减去等效节点力
[num,~] = size(DF);
for idf = 1:1:num
    if ie == DF( idf, 1 ) 
        enf = enf - BernoulliBeam2D_EquivalentNodeForce( DF(idf,1), ...
              DF(idf, 2), DF( idf, 3), DF( idf, 4 ) ) ;
        break ;
    end
end
% 把单元节点力转换到局部坐标下     
T = Beam2D_TransformMatrix( ie ) ;
enf = T  * enf ;
return