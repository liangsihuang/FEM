function enf = NonlinearBeam2D_ElementNodeForce( ie )
% 输入:
% ie:节点号
% 输出：
% enf：单元局部坐标系下的节点力 
global gElement gDelta
i = gElement( ie, 1 ) ;
j = gElement( ie, 2 ) ;
% 提取单元位移，注意：位移是整体坐标！
d_element = zeros( 6, 1 ) ;
d_element( 1:3 ) = gDelta( (i-1)*3+1:(i-1)*3+3 ) ;
d_element( 4:6 ) = gDelta( (j-1)*3+1:(j-1)*3+3 ) ;
k = NonlinearBeam2D_Stiffness(ie,1);%1代表整体坐标，因为位移是整体坐标的
enf = k * d_element;
% 先假设没有等效节点力......

% 把单元节点力转换到局部坐标下     
T = Beam2D_TransformMatrix( ie ) ;
enf = transpose( T ) * enf ;
return