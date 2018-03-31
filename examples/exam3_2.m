function exam3_2
% 本程序为第三章的第二个算例，采用平面梁单元计算斜腿刚架桥的变形和内力
%   输入参数： 无
%   输出结果： 节点位移和单元节点力 

    % 检查数据文件是否存在
    file_in = 'exam3_2.dat' ;
    if exist( file_in ) == 0
        disp( sprintf( '错误：文件 %s 不存在', file_in ) )
        disp( sprintf( '程序终止' ) )
        return ;
    end

    % 读入模型，求解并显示结果
    PlaneFrameModel(file_in) ;   % 定义有限元模型
    SolveModel ;                 % 求解有限元模型
    DisplayResults ;             % 显示计算结果
return ;

function PlaneFrameModel( file_in )
%  定义平面杆系的有限元模型
%  输入参数：
%      file_in ------- 有限元模型数据文件
%  返回值：
%      无
%  说明：
%      该函数定义平面杆系的有限元模型数据：
%        gNode ------- 节点定义
%        gElement ---- 单元定义
%        gMaterial --- 材料定义，包括弹性模量，梁的截面积和梁的抗弯惯性矩
%        gBC1 -------- 第一类边界条件
%        gBC2 -------- 第二类约束条件, 处理铰接的结点

    global gNode gElement gMaterial gBC1 gBC2

    % 打开数据文件
    fid_in = fopen( file_in, 'r' ) ;
    
    % 读取节点坐标和节点温度
    node_number = fscanf( fid_in, '%d', 1 ) ;
    gNode = zeros( node_number, 3 ) ;
    for i=1:node_number
        dummy = fscanf( fid_in, '%d', 1 ) ;
        gNode(i,:) = fscanf( fid_in, '%f', [1 3] ) ;
    end
     
    % 读取单元定义
    element_number = fscanf( fid_in, '%d', 1 ) ;
    gElement = zeros( element_number, 3 ) ;
    for i=1:element_number
        dummy = fscanf( fid_in, '%d', 1 ) ;
        gElement(i,:) = fscanf( fid_in, '%d', [1 3] ) ;
    end
        
    % 读取材料性质 
    material_number = fscanf( fid_in, '%d', 1 ) ;
    gMaterial = zeros( material_number, 5 ) ;
    for i=1:material_number
        dummy = fscanf( fid_in, '%d', 1 ) ;
        gMaterial(i,:) = fscanf( fid_in, '%f', [1 5] ) ;
    end
    
    % 读取第一类约束条件
    bc1_number = fscanf( fid_in, '%d', 1 ) ;
    gBC1 = zeros( bc1_number, 3 ) ;
    for i=1:bc1_number
        gBC1(i,1) = fscanf( fid_in, '%d', 1 ) ;
        gBC1(i,2) = fscanf( fid_in, '%d', 1 ) ;
        gBC1(i,3) = fscanf( fid_in, '%f', 1 ) ;
    end

    % 读取第二类约束条件
    bc2_number = fscanf( fid_in, '%d', 1 ) ;
    gBC2 = zeros( bc2_number, 3 ) ;
    for i=1:bc2_number
        gBC2(i,:) = fscanf( fid_in, '%d', [1 3] ) ;
    end
    
    % 关闭数据文件
    fclose( fid_in ) ;
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

    global gNode gElement gMaterial gBC1 gBC2 gK gDelta

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

    % step3. 计算自重温度变化产生的等效节点力
    for ie=1:1:element_number
        egf = EquivalentGravityForce( ie ) ;
        etf = EquivalentThermalForce( ie ) ;
        i = gElement( ie, 1 ) ;
        j = gElement( ie, 2 ) ;
        f( (i-1)*3+1 : (i-1)*3+3 ) = f( (i-1)*3+1 : (i-1)*3+3 ) + egf( 1:3 ) + etf( 1:3 );
        f( (j-1)*3+1 : (j-1)*3+3 ) = f( (j-1)*3+1 : (j-1)*3+3 ) + egf( 4:6 ) + etf( 4:6 );
    end
    
    % step4. 处理第二类约束条件，修改刚度矩阵和节点力向量。（采用划行划列法）
    [bc2_number,dummy] = size( gBC2 ) ;
    for ibc=1:1:bc2_number
        n1 = gBC2(ibc, 1 ) ;
        n2 = gBC2(ibc, 2 ) ;
        d = gBC2(ibc, 3 ) ;
        m1 = (n1-1)*3 + d ;
        m2 = (n2-1)*3 + d ;
        
        f(m2) = f(m2)+f(m1) ;
        f(m1) = 0 ;
        
        gK(m2,:) = gK(m2,:)+gK(m1,:) ;
        gK(:,m2) = gK(:,m2)+gK(:,m1) ;
        
        gK(m1,:) = zeros(1, node_number*3) ;
        gK(:,m1) = zeros(node_number*3, 1) ;
        gK(m1,m1) = 1.0 ;
        gK(m1,m2) = -1 ;
        gK(m2,m1) = -1 ;
        gK(m2,m2) = gK(m2,m2)+1 ;
    end
    
    % step5. 处理第一类约束条件，修改刚度矩阵和节点力向量。（采用划行划列法）
    [bc1_number,dummy] = size( gBC1 ) ;
    for ibc=1:1:bc1_number
        n = gBC1(ibc, 1 ) ;
        d = gBC1(ibc, 2 ) ;
        m = (n-1)*3 + d ;
        f = f - gBC1(ibc,3) * gK(:,m) ;
        f(m) = gBC1(ibc,3) ;
        gK(:,m) = zeros( node_number*3, 1 ) ;
        gK(m,:) = zeros( 1, node_number*3 ) ;
        gK(m,m) = 1.0 ;
    end

    % step 7. 求解方程组，得到节点位移向量
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

function egf = EquivalentGravityForce( ie )
%  计算单元自重的等效节点力
%  输入参数
%      ie  ----- 节点号
%  返回值
%      egf ----- 整体坐标系下的等效节点力 
    global gElement gNode gMaterial
    
    xi = gNode( gElement( ie, 1 ), 1 ) ;
    yi = gNode( gElement( ie, 1 ), 2 ) ;
    xj = gNode( gElement( ie, 2 ), 1 ) ;
    yj = gNode( gElement( ie, 2 ), 2 ) ;
    L = sqrt( (xj-xi)^2 + (yj-yi)^2 ) ;
    c = (xj-xi)/L ;
    
    A  = gMaterial( gElement( ie, 3 ), 3 ) ;
    ro = gMaterial( gElement( ie, 3 ), 4 ) ;
    g = 9.8 ;
    
    egf = [ 0; -ro*g*A*L/2; -ro*g*A*L^2*c/12; 0; -ro*g*A*L/2; ro*g*A*L^2*c/12 ] ;
return

function etf = EquivalentThermalForce( ie )
%  计算单元的温度荷载的等效节点力
%  输入参数
%      ie  ----- 节点号
%  返回值
%      etf ----- 整体坐标系下的等效节点力 
    global gElement gNode gMaterial
    
    dT1 = gNode( gElement( ie, 1 ), 3 ) ;
    dT2 = gNode( gElement( ie, 2 ), 3 ) ;
    E = gMaterial( gElement( ie, 3 ), 1 ) ;
    A = gMaterial( gElement( ie, 3 ), 3 ) ;
    alpha = gMaterial( gElement( ie, 3 ), 5 ) ;
    Nx = E*A*alpha ;
    
    etf = [ -Nx; 0; 0; Nx; 0; 0 ] ;
    T = TransformMatrix( ie ) ;
    etf = T * etf ;
    etf = (dT1+dT2)/2 * etf ;
return

function enf = ElementNodeForce( ie )
%  计算单元的节点力
%  输入参数
%      ie  ----- 节点号
%  返回值
%      enf ----- 单元局部坐标系下的节点力 
    global gElement gNode gDelta
    i = gElement( ie, 1 ) ;
    j = gElement( ie, 2 ) ;
    de = zeros( 6, 1 ) ;
    de( 1:3 ) = gDelta( (i-1)*3+1:(i-1)*3+3 ) ;
    de( 4:6 ) = gDelta( (j-1)*3+1:(j-1)*3+3 ) ;
    k = StiffnessMatrix( ie, 1 ) ;
    enf = k * de ;
    etf = EquivalentThermalForce( ie ) ;
    egf = EquivalentGravityForce( ie ) ;
    enf = enf - etf -egf ;
    T = TransformMatrix( ie ) ;
    enf = transpose( T ) * enf ;
return