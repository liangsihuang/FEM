function BernoulliBeam2D_SolveModel
%      该函数求解有限元模型，过程如下
%        1. 计算单元刚度矩阵，集成整体刚度矩阵
%        2. 计算单元的等效节点力，集成整体节点力向量
%        3. 处理约束条件，修改整体刚度矩阵和节点力向量
%        4. 求解方程组，得到整体节点位移向量

    global gNode gElement gMaterial gBC1 gNF gDF gK gDelta

    % step1. 定义整体刚度矩阵和节点力向量
    [node_number,~] = size( gNode ) ;
%     sparse 定义稀疏矩阵,gDelta = gK \ f也是稀疏矩阵
    gK = sparse( node_number * 3, node_number * 3 ) ;
    f = sparse( node_number * 3, 1 ) ;

    % step2. 计算单元刚度矩阵，并集成到整体刚度矩阵中
    [element_number,~] = size( gElement ) ;
    for ie=1:1:element_number
        k = BernoulliBeam2D_Stiffness( ie, 1 ) ;
        BernoulliBeam2D_Assembly( ie, k ) ;
    end

    % step3. 把集中力直接集成到整体节点力向量中
    [nf_number, ~] = size( gNF ) ;
    for inf=1:1:nf_number
        n = gNF( inf, 1 ) ;
        d = gNF( inf, 2 ) ;
        f( (n-1)*3 + d ) = gNF( inf, 3 ) ;
    end

    % step4. 计算分布力的等效节点力，集成到整体节点力向量中
    [df_number, ~] = size( gDF ) ;
    for idf = 1:1:df_number
        enf = BernoulliBeam2D_EquivalentNodeForce( gDF(idf,1), gDF(idf, 2), gDF( idf, 3), gDF( idf, 4 ) ) ;
        i = gElement( gDF(idf,1), 1 ) ;
        j = gElement( gDF(idf,1), 2 ) ;
        f( (i-1)*3+1 : (i-1)*3+3 ) = f( (i-1)*3+1 : (i-1)*3+3 ) + enf( 1:3 ) ;
        f( (j-1)*3+1 : (j-1)*3+3 ) = f( (j-1)*3+1 : (j-1)*3+3 ) + enf( 4:6 ) ;
    end
  
    % step5. 处理约束条件，修改刚度矩阵和节点力向量。采用乘大数法
    [bc_number,~] = size( gBC1 ) ;
    for ibc=1:1:bc_number
        n = gBC1(ibc, 1 ) ;
        d = gBC1(ibc, 2 ) ;
        m = (n-1)*3 + d ;
        f(m) = gBC1(ibc, 3)* gK(m,m) * 1e15 ;
        gK(m,m) = gK(m,m) * 1e15 ;
    end

    % step 6. 求解方程组，得到节点位移向量
    gDelta = gK \ f ;
return

function k=BernoulliBeam2D_Stiffness(ie,icoord)
%  计算单元刚度矩阵
%  输入参数:
%     ie -------  单元号
%     icoord  --  坐标系参数，可以是下面两个之一
%                    1  ----  整体坐标系
%                    2  ----  局部坐标系
%  返回值:
%     k  ----  根据icoord的值，相应坐标系下的刚度矩阵
    global gNode gElement gMaterial
    k = zeros( 6, 6 ) ;
    E = gMaterial( gElement(ie, 3), 1 ) ;
    I = gMaterial( gElement(ie, 3), 2 ) ;
    A = gMaterial( gElement(ie, 3), 3 ) ;
    xi = gNode( gElement( ie, 1 ), 1 ) ;
    yi = gNode( gElement( ie, 1 ), 2 ) ;
    xj = gNode( gElement( ie, 2 ), 1 ) ;
    yj = gNode( gElement( ie, 2 ), 2 ) ;
    L = ( (xj-xi)^2 + (yj-yi)^2 )^(1/2) ;
    k = [  E*A/L           0          0 -E*A/L           0          0
               0  12*E*I/L^3  6*E*I/L^2      0 -12*E*I/L^3  6*E*I/L^2
               0   6*E*I/L^2    4*E*I/L      0  -6*E*I/L^2    2*E*I/L
          -E*A/L           0          0  E*A/L           0          0
               0 -12*E*I/L^3 -6*E*I/L^2      0  12*E*I/L^3 -6*E*I/L^2
               0   6*E*I/L^2    2*E*I/L      0  -6*E*I/L^2    4*E*I/L] ;
    if icoord == 1
        T = BernoulliBeam2D_TransformMatrix( ie ) ;
        k = T*k*transpose(T) ;
    end
return

function BernoulliBeam2D_Assembly( ie, k )
%  把单元刚度矩阵集成到整体刚度矩阵
%  输入参数:
%      ie  --- 单元号
%      k   --- 单元刚度矩阵
%  返回值:
%      无
    global gElement gK
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

function T = BernoulliBeam2D_TransformMatrix( ie )
%  计算单元的坐标转换矩阵( 局部坐标 -> 整体坐标 )
%  输入参数
%      ie  ----- 节点号
%  返回值
%      T ------- 从局部坐标到整体坐标的坐标转换矩阵
    global gElement gNode
    xi = gNode( gElement( ie, 1 ), 1 ) ;
    yi = gNode( gElement( ie, 1 ), 2 ) ;
    xj = gNode( gElement( ie, 2 ), 1 ) ;
    yj = gNode( gElement( ie, 2 ), 2 ) ;
    L = sqrt( (xj-xi)^2 + (yj-yi)^2 ) ;
    c = (xj-xi)/L ;
    s = (yj-yi)/L ;
    T=[ c  -s   0   0   0   0
        s   c   0   0   0   0
        0   0   1   0   0   0
        0   0   0   c  -s   0
        0   0   0   s   c   0
        0   0   0   0   0   1] ;
return

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