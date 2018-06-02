function T =Beam2D_TransformMatrix( ie )
%  计算单元的坐标转换矩阵( 局部坐标 -> 整体坐标 )
%  输入参数
%      ie  ----- 节点号
%  返回值
%      T ------- 从局部坐标到整体坐标的坐标转换矩阵
global Element Node
xi = Node( Element( ie, 1 ), 1 ) ;
yi = Node( Element( ie, 1 ), 2 ) ;
xj = Node( Element( ie, 2 ), 1 ) ;
yj = Node( Element( ie, 2 ), 2 ) ;
L = sqrt( (xj-xi)^2 + (yj-yi)^2 ) ;
c = (xj-xi)/L ;
s = (yj-yi)/L ;
T=[ c  s   0   0   0   0
    -s   c   0   0   0   0
    0   0   1   0   0   0
    0   0   0   c  s   0
    0   0   0   -s   c   0
    0   0   0   0   0   1] ;
return