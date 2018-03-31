function exam3_1
% 本程序为第三章的第一个算例，采用平面梁单元计算平面刚架的变形和内力
%      输入参数： 无
%      输出结果： 节点位移和单元节点力 

    PlaneFrameModel ;      % 定义有限元模型
    SolveModel ;           % 求解有限元模型
    DisplayResults ;       % 显示计算结果
return ;

function PlaneFrameModel
%  定义平面杆系的有限元模型
%  输入参数：
%      无
%  返回值：
%      无
%  说明：
%      该函数定义平面杆系的有限元模型数据：
%        gNode ------- 节点定义
%        gElement ---- 单元定义
%        gMaterial --- 材料定义，包括弹性模量，梁的截面积和梁的抗弯惯性矩
%        gBC1 -------- 约束条件
%        gNF --------- 集中力
%        gDF --------- 分布力

    global gNode gElement gMaterial gBC1 gNF gDF

    % 节点坐标
    %         x      y
    gNode = [0.0,   0.0          % 节点 1
             0.0,   4.0          % 节点 2
             3.0,   0.0          % 节点 3
             3.0,   4.0          % 节点 4
             4.5,   4.0          % 节点 5
             6.0,   0.0          % 节点 6
             6.0,   4.0 ] ;      % 节点 7
     
    % 单元定义
    %          节点1  节点2  材料号
    gElement = [1,      2,      1       % 单元 1
                2,      4,      1       % 单元 2
                3,      4,      1       % 单元 3
                4,      5,      1       % 单元 4
                5,      7,      1       % 单元 5
                6,      7,      1] ;    % 单元 6
        
    % 材料性质 
    %           弹性模量  抗弯惯性矩  截面积
    gMaterial = [2.1e11,   2.0e-4,    1.0e-2] ;   %  材料 1

    % 第一类约束条件
    %     节点号   自由度号    约束值
    gBC1 = [ 1,        1,        0.0
             1,        2,        0.0
             1,        3,        0.0
             3,        1,        0.0
             3,        2,        0.0
             3,        3,        0.0
             6,        1,        0.0
             6,        2,        0.0
             6,        3,        0.0] ;

    % 集中力
    %     节点号   自由度号   集中力值
    gNF = [  5,       2,         -80e3] ;

    % 分布载荷（线性分布）
    %     单元号   节点1载荷值   节点2载荷值   自由度号
    gDF = [  1         -30e3         0              2
             2         -15e3       -15e3            2  ] ;
return

function SolveModel
%  求解有限元模型
%  输入参数：
%     无
%  返回值：
%     无
%  说明：
%      该函数求解有限元模型，过程如下
%        1. 计算单元刚度矩阵，集成整体刚度矩阵
%        2. 计算单元的等效节点力，集成整体节点力向量
%        3. 处理约束条件，修改整体刚度矩阵和节点力向量
%        4. 求解方程组，得到整体节点位移向量

    global gNode gElement gMaterial gBC1 gNF gDF gK gDelta

    % step1. 定义整体刚度矩阵和节点力向量
    [node_number,dummy] = size( gNode ) ;
    gK = sparse( node_number * 3, node_number * 3 ) ;
    f = sparse( node_number * 3, 1 ) ;

    % step2. 计算单元刚度矩阵，并集成到整体刚度矩阵中
    [element_number,dummy] = size( gElement ) ;
    for ie=1:1:element_number
        k = StiffnessMatrix( ie, 1 ) ;
        AssembleStiffnessMatrix( ie, k ) ;
    end

    % step3. 把集中力直接集成到整体节点力向量中
    [nf_number, dummy] = size( gNF ) ;
    for inf=1:1:nf_number
        n = gNF( inf, 1 ) ;
        d = gNF( inf, 2 ) ;
        f( (n-1)*3 + d ) = gNF( inf, 3 ) ;
    end

    % step4. 计算分布力的等效节点力，不集成到整体节点力向量中
    [df_number, dummy] = size( gDF ) ;
    for idf = 1:1:df_number
        enf = EquivalentNodeForce( gDF(idf,1), gDF(idf, 2), gDF( idf, 3), gDF( idf, 4 ) ) ;
        i = gElement( gDF(idf,1), 1 ) ;
        j = gElement( gDF(idf,1), 2 ) ;
        f( (i-1)*3+1 : (i-1)*3+3 ) = f( (i-1)*3+1 : (i-1)*3+3 ) + enf( 1:3 ) ;
        f( (j-1)*3+1 : (j-1)*3+3 ) = f( (j-1)*3+1 : (j-1)*3+3 ) + enf( 4:6 ) ;
    end
  
    % step5. 处理约束条件，修改刚度矩阵和节点力向量。采用乘大数法
    [bc_number,dummy] = size( gBC1 ) ;
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

function DisplayResults
%  显示计算结果
%  输入参数：
%     无
%  返回值：
%     无

    global gNode gElement gMaterial gBC1 gNF gDF gK gDelta
    
    fprintf( '节点位移\n' ) ; 
    fprintf( '  节点号         x方向位移               y方向位移               转角\n' ) ; 
    [node_number,dummy] = size( gNode ) ;
    for i=1:node_number
        fprintf(  '%6d       %16.8e        %16.8e       %16.8e\n',...
                  i, gDelta((i-1)*3+1), gDelta((i-1)*3+2), gDelta((i-1)*3+3) ) ; 
    end
    fprintf( '\n\n节点力\n' ) ; 
    fprintf( '                                        轴力               剪力               弯矩\n' ) ; 
    [element_number, dummy] = size( gElement ) ;
    for ie = 1:element_number
        enf = ElementNodeForce( ie ) ;
        fprintf( '单元号%6d    节点号%6d     %16.8e  %16.8e  %16.8e\n', ...
                  ie, gElement(ie,1), enf(1), enf(2), enf(3) ) ;
        fprintf( '                节点号%6d     %16.8e  %16.8e  %16.8e\n', ...
                  gElement(ie,2), enf(4), enf(5), enf(6) ) ;
    end
return

function k = StiffnessMatrix( ie, icoord )
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
        T = TransformMatrix( ie ) ;
        k = T*k*transpose(T) ;
    end
return

function AssembleStiffnessMatrix( ie, k )
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

function enf = EquivalentNodeForce( ie, p1, p2, idof )
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
    enf = zeros( 6, 1 ) ;                       % 定义 6x1 的等效节点力向量
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
    
    T = TransformMatrix( ie ) ;             %  计算单元的转换矩阵
    enf = T * enf ;                         %  把等效节点力转换到整体坐标下
return

function T = TransformMatrix( ie )
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

function enf = ElementNodeForce( ie )
%  计算单元的节点力
%  输入参数
%      ie  ----- 节点号
%  返回值
%      enf ----- 单元局部坐标系下的节点力 
    global gElement gNode gDelta gDF
    i = gElement( ie, 1 ) ;
    j = gElement( ie, 2 ) ;
    de = zeros( 6, 1 ) ;
    de( 1:3 ) = gDelta( (i-1)*3+1:(i-1)*3+3 ) ;
    de( 4:6 ) = gDelta( (j-1)*3+1:(j-1)*3+3 ) ;
    k = StiffnessMatrix( ie, 1 ) ;
    enf = k * de ;
    
    [df_number, dummy] = size( gDF ) ;
    for idf = 1:1:df_number
        if ie == gDF( idf, 1 ) 
            enf = enf - EquivalentNodeForce( gDF(idf,1), ...
                  gDF(idf, 2), gDF( idf, 3), gDF( idf, 4 ) ) ;
            break ;
        end
    end
    
    T = TransformMatrix( ie ) ;
    enf = transpose( T ) * enf ;
return