function exam4_1( m, n )
% 本程序为第四章的第一个算例，采用矩形单元计算纯弯梁的变形
%      exam4_1(m,n) 
%  输入参数： 
%      m ------  x方向单元数目
%      n ------  y方向单元数目

% 定义全局变量
%      gNode ------ 节点坐标
%      gElement --- 单元定义
%      gMaterial -- 材料性质
%      gBC1 ------- 第一类约束条件
%      gDF -------- 分布力
%      gK --------- 整体刚度矩阵
%      gDelta ----- 整体节点坐标
    global gNode gElement gMaterial gBC1 gDF gK gDelta

    if nargin < 1 
        m = 4 ;
        n = 4 ;
    elseif nargin < 2
        n = 4 ;
    end
    
    FemModel(m, n) ;       % 定义有限元模型
    SolveModel ;           % 求解有限元模型
    DisplayResults ;       % 显示计算结果
return ;

function FemModel(m, n)
%  定义有限元模型
%  输入参数：
%      m ---  x方向单元数目
%      n ---  y方向单元数目
%  返回值：
%      无
%  说明：
%      该函数定义平面杆系的有限元模型数据：
%        gNode ------- 节点定义
%        gElement ---- 单元定义
%        gMaterial --- 材料定义，包括弹性模量，梁的截面积和梁的抗弯惯性矩
%        gBC --------- 约束条件
%        gDF --------- 分布力

    global gNode gElement gMaterial gBC1 gDF

    length = 0.045 ;    % 计算部分的长度(x方向)
    height = 0.03 ;     % 计算部分的高度(y方向)
    dx = length/m ;     % 矩形单元的宽度
    dy = height/n ;     % 矩形单元的高度
    p = 10e6 ;          % 弯曲应力10MPa
    
    % 节点坐标
    gNode = zeros( (m+1)*(n+1), 2 ) ;
    for i=1:n+1
        for j=1:m+1
            k = (i-1)*(m+1)+j ;        % 节点号
            xk = (j-1)*dx ;            % 节点的x坐标
            yk = (i-1)*dy ;            % 节点的y坐标
            gNode(k,:) = [xk, yk] ;
        end
    end
     
    % 单元定义
    gElement = zeros( m*n, 4 ) ;
    for i=1:n
        for j=1:m
            k = (i-1)*m+j ;             % 单元号
            n1 = (i-1)*(m+1)+j ;        % 第一个节点号
            n2 = (i-1)*(m+1)+j+1 ;      % 第二个节点号
            n3 = i*(m+1)+j+1 ;          % 第三个节点号
            n4 = i*(m+1)+j ;            % 第四个节点号
            gElement(k,:) = [n1, n2, n3, n4] ;
        end
    end

    % 材料性质 
    %           弹性模量    泊松比   厚度
    gMaterial = [2.0e11,    0.3,     0.01] ;   %  材料 1

    % 第一类约束条件
    gBC1 = zeros( m+1+n+1, 3 ) ;
    for j=1:m+1
        gBC1(j,:) = [j, 1, 0.0] ;
    end
    for i=2:n+1
        gBC1(m+1+i,:) = [(i-1)*(m+1)+1, 1, 0.0] ; 
    end
    gBC1(m+1+1,:) = [1, 2, 0.0] ;

    % 分布载荷（线性分布）
    gDF = zeros( n, 5 ) ;
    for i=1:n
        k = i*m ;
        gDF(i,:) = [ k, 2, (i-1)*p/n, i*p/n, 1] ;
    end
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

    global gNode gElement gMaterial gBC1 gDF gK gDelta

    % step1. 定义整体刚度矩阵和节点力向量
    [node_number,dummy] = size( gNode ) ;
    gK = sparse( node_number * 2, node_number * 2 ) ;
    f = sparse( node_number * 2, 1 ) ;

    % step2. 计算单元刚度矩阵，并集成到整体刚度矩阵中
    [element_number,dummy] = size( gElement ) ;
    for ie=1:1:element_number
        k = StiffnessMatrix( ie ) ;
        AssembleStiffnessMatrix( ie, k ) ;
    end

    % step3. 计算分布力的等效节点力，并集成到整体节点力向量中
    [df_number, dummy] = size( gDF ) ;
    for idf = 1:1:df_number
        enf = EquivalentNodeForce( gDF(idf,1), gDF(idf,2), gDF(idf,3), gDF(idf,4), gDF(idf,5) ) ;
        ielem = gDF(idf,1) ;
        iedge = gDF(idf,2) ;
        i = gElement( ielem, iedge ) ;
        if iedge < 4 
            j = gElement( ielem, iedge+1 ) ;
        else
            j = gElement( ielem, 1 );
        end
        
        f( (i-1)*2+1 : (i-1)*2+2 ) = f( (i-1)*2+1 : (i-1)*2+2 ) + enf( 1:2 ) ;
        f( (j-1)*2+1 : (j-1)*2+2 ) = f( (j-1)*2+1 : (j-1)*2+2 ) + enf( 3:4 ) ;
    end
  
    % step4. 处理约束条件，修改刚度矩阵和节点力向量。采用乘大数法
    [bc_number,dummy] = size( gBC1 ) ;
    for ibc=1:1:bc_number
        n = gBC1(ibc, 1 ) ;
        d = gBC1(ibc, 2 ) ;
        m = (n-1)*2 + d ;
        f(m) = gBC1(ibc, 3)* gK(m,m) * 1e15 ;
        gK(m,m) = gK(m,m) * 1e15 ;
    end

    % step 5. 求解方程组，得到节点位移向量
    gDelta = gK \ f ;
return


function k = StiffnessMatrix( ie )
%  计算单元刚度矩阵
%  输入参数:
%     ie ----  单元号
%  返回值:
%     k  ----  单元刚度矩阵

global gNode gElement gMaterial 
    k = zeros( 8, 8 ) ;
    E  = gMaterial( 1 ) ;
    mu = gMaterial( 2 ) ;
    h  = gMaterial( 3 ) ;
    x1 = gNode( gElement( ie, 1 ), 1 ) ;
    y1 = gNode( gElement( ie, 1 ), 2 ) ;
    x3 = gNode( gElement( ie, 3 ), 1 ) ;
    y3 = gNode( gElement( ie, 3 ), 2 ) ;
    a = (x3-x1)/2 ;
    b = (y3-y1)/2 ;
    xi  = [ -1, 1, 1, -1 ] ;
    eta = [ -1, -1, 1, 1 ] ;
    for i=1:1:4
        for j=1:1:4
            k( (i-1)*2 + 1, (j-1)*2 + 1 ) = b/a*xi(i)*xi(j)*(1+1/3*eta(i)*eta(j)) + ...
                (1-mu)/2*a/b*eta(i)*eta(j)*(1+1/3*xi(i)*xi(j)); 
            k( (i-1)*2 + 1, (j-1)*2 + 2 ) = mu*eta(j)*xi(i) + (1-mu)/2*xi(j)*eta(i) ;
            k( (i-1)*2 + 2, (j-1)*2 + 1 ) = mu*xi(j)*eta(i) + (1-mu)/2*eta(j)*xi(i) ;
            k( (i-1)*2 + 2, (j-1)*2 + 2 ) = a/b*eta(i)*eta(j)*(1+1/3*xi(i)*xi(j)) + ...
                (1-mu)/2*b/a*xi(i)*xi(j)*(1+1/3*eta(i)*eta(j)) ;
        end
    end
    eh = E*h/4/(1-mu^2) ;
    k = eh*k ;
return

function AssembleStiffnessMatrix( ie, k )
%  把单元刚度矩阵集成到整体刚度矩阵
%  输入参数:
%      ie  --- 单元号
%      k   --- 单元刚度矩阵
%  返回值:
%      无
    global gElement gK
    for i=1:1:4
        for j=1:1:4
            for p=1:1:2
                for q=1:1:2
                    m = (i-1)*2+p ;
                    n = (j-1)*2+q ;
                    M = (gElement(ie,i)-1)*2+p ;
                    N = (gElement(ie,j)-1)*2+q ;
                    gK(M,N) = gK(M,N) + k(m,n) ;
                end
            end
        end
    end
return

function enf = EquivalentNodeForce( ie, iedge, p1, p2, idof )
%   计算线性分布荷载的等效节点力
%   输入参数:
%      ie  -----  单元号
%      ieged ---  作用的边号
%      p1  -----  第一个节点上的分布力集度值
%      p2  -----  第二个节点上的分布力集度值
%      idof  ---  分布力的方向
%                  1 ---  x方向
%                  2 ---  y方向
%   返回值:
%      enf -----  等效节点力向量
    global gElement gNode gMaterial
    h = gMaterial( 3 ) ;
    x1 = gNode( gElement( ie, 1 ), 1 ) ;
    y1 = gNode( gElement( ie, 1 ), 2 ) ;
    x3 = gNode( gElement( ie, 3 ), 1 ) ;
    y3 = gNode( gElement( ie, 3 ), 2 ) ;
    a = ( x3 - x1 ) / 2 ;
    b = ( y3 - y1 ) / 2 ;
    if iedge == 1 | iedge == 3
        length = a ;        
    else
        length = b ;
    end
    f1 = h*length*(2*p1+p2)/3 ;
    f2 = h*length*(p1+2*p2)/3 ;
    if idof == 1
        enf = [f1; 0; f2; 0] ;
    else
        enf = [0; f1; 0; f2] ;
    end
return


function DisplayResults
%  显示计算结果
%  输入参数：
%     无
%  返回值：
%     无

    global gNode gDelta
    
    fprintf( '节点位移\n' ) ; 
    fprintf( '  节点号         x方向位移               y方向位移\n' ) ; 
    [node_number,dummy] = size( gNode ) ;
    for i=1:node_number
        fprintf(  '%6d       %16.8e        %16.8e\n',...
                  i, gDelta((i-1)*2+1), gDelta((i-1)*2+2) ) ; 
    end
return
