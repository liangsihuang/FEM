function enf = BernoulliBeam2D_EquivalentNodeForce( ie, p1, p2, idof )
%   计算线性分布荷载的等效节点力
%   输入参数:
%      ie  -----  单元号
%      p1  -----  第一个节点上的分布力集度值
%      p2  -----  第二个节点上的分布力集度值
%      idof  ---  分布力的种类，它可以是下面几种
%                  1 ---  分布轴向力
%                  2 ---  分布横向力
%                  3 ---  分布弯矩
%   返回值:
%      enf -----  整体坐标系下等效节点力向量
    global gElement gNode
    enf = zeros( 6, 1 ) ; % 定义 6x1 的等效节点力向量
    xi = gNode( gElement( ie, 1 ), 1 ) ;
    yi = gNode( gElement( ie, 1 ), 2 ) ;
    xj = gNode( gElement( ie, 2 ), 1 ) ;
    yj = gNode( gElement( ie, 2 ), 2 ) ;
    L = sqrt( (xj-xi)^2 + (yj-yi)^2 ) ;
    switch idof 
    case 1     %  分布轴向力 
        enf( 1 ) = (2*p1+p2)*L/6 ;
        enf( 4 ) = (p1+2*p2)*L/6 ;
    case 2     %  分布横向力
        enf( 2 ) = (7*p1+3*p2)*L/20 ;
        enf( 3 ) = (3*p1+2*p2)*L^2/60 ;
        enf( 5 ) = (3*p1+7*p2)*L/20 ;
        enf( 6 ) = -(2*p1+3*p2)*L^2/60 ;
    case 3     %  分布弯矩
        enf( 2 ) = -(p1+p2)/2 ;
        enf( 3 ) = (p1-p2)*L/12 ;
        enf( 5 ) = (p1+p2)/2 ;
        enf( 6 ) = -(p1-p2)*L/12 ;
    otherwise
        disp( sprintf( '分布力的种类错误，单元号:%d',ie ) ) ;
    end
    
    T = BernoulliBeam2D_TransformMatrix( ie );%  计算单元的转换矩阵
    enf = T * enf;%  把等效节点力转换到整体坐标下
return