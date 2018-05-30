function enf = BernoulliBeam2D_EquivalentNodeForce( ie, p1, p2, force_type )
% 输入:
% ie:单元号
% p1:第一个节点上的分布力集度值
% p2:第二个节点上的分布力集度值
% force_type:分布力的种类:1代表分布轴向力，2代表分布横向力，3代表分布弯矩
% 输出:
% enf:整体坐标系下等效节点力向量
global Element Node
enf = zeros( 6, 1 ) ; % 定义 6x1 的等效节点力向量
xi = Node( Element( ie, 1 ), 1 ) ;
yi = Node( Element( ie, 1 ), 2 ) ;
xj = Node( Element( ie, 2 ), 1 ) ;
yj = Node( Element( ie, 2 ), 2 ) ;
L = sqrt( (xj-xi)^2 + (yj-yi)^2 ) ;
if force_type==1 
    % 分布轴向力 
    enf( 1 ) = (2*p1+p2)*L/6 ;
    enf( 4 ) = (p1+2*p2)*L/6 ;
end
if force_type==2
    % 分布横向力
    enf( 2 ) = (7*p1+3*p2)*L/20 ;
    enf( 3 ) = (3*p1+2*p2)*L^2/60 ;
    enf( 5 ) = (3*p1+7*p2)*L/20 ;
    enf( 6 ) = -(2*p1+3*p2)*L^2/60 ;
end
if force_type==3
    % 分布弯矩
    enf( 2 ) = -(p1+p2)/2 ;
    enf( 3 ) = (p1-p2)*L/12 ;
    enf( 5 ) = (p1+p2)/2 ;
    enf( 6 ) = -(p1-p2)*L/12 ;
end
%  把等效节点力转换到整体坐标下
T = Beam2D_TransformMatrix( ie );
enf = T' * enf;
return