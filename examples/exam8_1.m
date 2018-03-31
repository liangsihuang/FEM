function exam8_1
% 本程序为第八章的第一个算例，采用平面梁单元计算两铰抛物线拱的自由振动特性
%      输入参数： 无
%      输出结果： 前3阶振动频率及其相应的振型 

% 定义全局变量
%      gNode ------ 节点坐标
%      gElement --- 单元定义
%      gMaterial -- 材料性质
%      gBC1 ------- 第一类约束条件
%      gK --------- 整体刚度矩阵
%      gDelta ----- 整体节点坐标

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
%        gBC --------- 约束条件

    global gNode gElement gMaterial gBC1

    % 给定抛物线拱的几何特征
    L = 60 ;               %  计算跨径     
    f = 7.5 ;              %  计算矢高
    
    n = 100 ;              %  单元数目
    x = -L/2:L/n:L/2 ;     %  结点的x坐标
    a = f/L^2*4 ;
    y = - a * x.^2 ;       %  结点的y坐标

    % 节点坐标
    gNode = [x'  y'] ;
    
    % 单元定义
    gElement = zeros( n, 3 ) ;
    for i=1:n
        gElement( i, : ) = [ i, i+1, 1 ] ;
    end
    
    % 材料性质 
    %           弹性模量   抗弯惯性矩   截面积   密度
    gMaterial = [2.06e11,  0.03622,   0.0815,  1435.2/0.0815];   %  材料 1

    % 第一类约束条件
    %     节点号   自由度号    约束值
    gBC1 = [ 1,        1,        0.0
             1,        2,        0.0
             n+1,      1,        0.0
             n+1,      2,        0.0] ;
return

function SolveModel
%  求解有限元模型
%  输入参数：
%     无
%  返回值：
%     无
%  说明：
%      该函数求解有限元模型，过程如下
%        1. 计算单元的刚度和质量矩阵，集成整体刚度和质量矩阵
%        2. 处理约束条件，修改整体刚度矩阵
%        3. 求解特征值问题

    global gNode gElement gMaterial gBC1 gK gM gEigValue gEigVector

    % step1. 定义整体刚度矩阵和节点力向量
    [node_number,dummy] = size( gNode ) ;
    gK = sparse( node_number * 3, node_number * 3 ) ;
    gM = sparse( node_number * 3, node_number * 3 ) ;

    % step2. 计算单元刚度和质量矩阵，并集成到整体刚度和质量矩阵中
    [element_number,dummy] = size( gElement ) ;
    for ie=1:1:element_number
        k = StiffnessMatrix( ie ) ;
        m = MassMatrix( ie ) ; 
        AssembleGlobalMatrix( ie, k, m ) ;
    end

    % step3. 处理第一类约束条件，修改刚度矩阵和质量矩阵。（采用划行划列法）
    [bc1_number,dummy] = size( gBC1 ) ;
    w2max = max( diag(gK)./diag(gM) ) ;
    for ibc=1:1:bc1_number
        n = gBC1(ibc, 1 ) ;
        d = gBC1(ibc, 2 ) ;
        m = (n-1)*3 + d ;
        gK(:,m) = zeros( node_number*3, 1 ) ;
        gK(m,:) = zeros( 1, node_number*3 ) ;
        gK(m,m) = 1;
        gM(:,m) = zeros( node_number*3, 1 ) ;
        gM(m,:) = zeros( 1, node_number*3 ) ;
        gM(m,m) = gK(m,m)/w2max/1e10 ;
    end
    
    % step4. 求解特征值问题
    % step4.1为了使刚度矩阵和质量矩阵对称（在计算时可能引入舍入误差）
    for i=1:node_number*3
        for j=i:node_number*3
            gK(j,i) = gK(i,j) ;
            gM(j,i) = gM(i,j) ;
        end
    end
    
    % step4.2 计算前6阶特征值和特征向量
    [gEigVector, gEigValue] = eigs(gK, gM, 3, 'SM' ) ; 
    
    % step4.3 修改特征向量中受约束的自由度
    for ibc=1:1:bc1_number
        n = gBC1(ibc, 1 ) ;
        d = gBC1(ibc, 2 ) ;
        m = (n-1)*3 + d ;
        gEigVector(m,:) = gBC1(ibc,3) ;
    end
return

function k = StiffnessMatrix( ie )
%  计算单元刚度矩阵
%  输入参数:
%     ie -------  单元号
%  返回值:
%     k  ----  整体坐标系下的刚度矩阵
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
    T = TransformMatrix( ie ) ;
    k = T*k*transpose(T) ;
return

function m = MassMatrix( ie )
%  计算单元质量矩阵
%  输入参数:
%     ie -------  单元号
%  返回值:
%     m  ----  整体坐标系下的质量矩阵
    global gNode gElement gMaterial
    m = zeros( 6, 6 ) ;
    E = gMaterial( gElement(ie, 3), 1 ) ;
    A = gMaterial( gElement(ie, 3), 3 ) ;
    ro = gMaterial( gElement(ie, 3 ), 4 ) ;
    xi = gNode( gElement( ie, 1 ), 1 ) ;
    yi = gNode( gElement( ie, 1 ), 2 ) ;
    xj = gNode( gElement( ie, 2 ), 1 ) ;
    yj = gNode( gElement( ie, 2 ), 2 ) ;
    L = ( (xj-xi)^2 + (yj-yi)^2 )^(1/2) ;
    m = ro*A*L/420*[140      0      0   70      0      0
                      0    156   22*L    0     54   -13*L
                      0   22*L  4*L^2    0   13*L  -3*L^2
                     70      0      0  140      0       0 
                      0     54   13*L    0    156   -22*L
                      0  -13*L -3*L^2    0  -22*L  4*L^2 ] ;
    T = TransformMatrix( ie ) ;
    m = T*m*transpose(T) ;
return

function AssembleGlobalMatrix( ie, ke, me )
%  把单元刚度和质量矩阵集成到整体刚度矩阵
%  输入参数:
%      ie  --- 单元号
%      ke  --- 单元刚度矩阵
%      me  --- 单元质量矩阵
%  返回值:
%      无
    global gElement gK gM
    for i=1:1:2
        for j=1:1:2
            for p=1:1:3
                for q =1:1:3
                    m = (i-1)*3+p ;
                    n = (j-1)*3+q ;
                    M = (gElement(ie,i)-1)*3+p ;
                    N = (gElement(ie,j)-1)*3+q ;
                    gK(M,N) = gK(M,N) + ke(m,n) ;
                    gM(M,N) = gM(M,N) + me(m,n) ;
                end
            end
        end
    end
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

function DisplayResults
%  显示计算结果
%  输入参数：
%     无
%  返回值：
%     无

    global gNode gElement gMaterial gBC1 gEigValue gEigVector

    fre_number = length(diag(gEigValue)) ;
    
    % 打印特征向量(振型)
    fprintf( '\n\n 表一   特征向量(振型)  \n' ) ;
    for i=1:fre_number
        fprintf( '----------------') ;
    end
    fprintf( '\n' ) ;
    for i=1:fre_number
        fprintf( '  %6d        ', i ) ;
    end
    fprintf( '\n' ) ;
    for i=1:fre_number
        fprintf( '----------------') ;
    end
    fprintf( '\n' ) ;
    [dof,dummy]=size(gEigVector) ;
    for i=1:dof
        for j=fre_number:-1:1
            fprintf( '%15.7e ', gEigVector(i,j) ) ;
        end
        fprintf( '\n' ) ;
    end
    for i=1:fre_number
        fprintf( '----------------') ;
    end
    fprintf( '\n' ) ;

    % 打印特征值
    fprintf( '\n\n\n\n 表二   特征值(频率)列表  \n' ) ;
    fprintf( '----------------------------------------------------------\n') ;
    fprintf( '   阶数        特征值          频率(Hz)        圆频率(Hz)  \n' ) ;
    fprintf( '----------------------------------------------------------\n') ;
    for i=fre_number:-1:1
        fprintf( '%6d   %15.7e   %15.7e   %15.7e\n', fre_number-i+1, ...
            gEigValue(i,i), sqrt(gEigValue(i,i))/2/pi, sqrt(gEigValue(i,i)) ) ;
    end
    fprintf( '----------------------------------------------------------\n') ;
    
    % 绘制振型图
    for j=fre_number:-1:1
        figure ;
        x = gNode(:,1) ;
        y = gNode(:,2) ;
        dx = gEigVector(1:3:length(x)*3, j ) ;
        dy = gEigVector(2:3:length(x)*3, j ) ;
        factor = max( [max(abs(x))/max(abs(dx)), max(abs(y))/max(abs(dy))] )* 0.05; 
        plot(x,y,'-', x+factor*dx, y+factor*dy, ':') ;
        title( sprintf( '第%d阶频率: %.3f Hz', fre_number-j+1, sqrt(gEigValue(j,j))/2/pi ) ) ;
        axis equal;
        axis off ;
    end
return
