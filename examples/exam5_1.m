function exam5_1
% 本程序为第五章的第一个算例，采用四面体单元计算单向拉伸杆的横截面应力
%  输入参数： 
%      无

% 定义全局变量
%      gNode ------ 节点坐标
%      gElement --- 单元定义
%      gMaterial -- 材料性质
%      gBC1 ------- 第一类约束条件
%      gDF -------- 分布力
%      gNF -------- 集中力
%      gK --------- 整体刚度矩阵
%      gDelta ----- 整体节点坐标

    file_in = 'exam5_1.dat' ;
    
    % 检查文件是否存在
    if exist( file_in ) == 0
        disp( sprintf( '错误：文件 %s 不存在', file_in ) )
        disp( sprintf( '程序终止' ) )
        return ;
    end
    
    % 根据输入文件名称生成输出文件名称
    [path_str,name_str,ext_str] = fileparts( file_in ) ;
    ext_str_out = '.mat' ;
    file_out = fullfile( path_str, [name_str, ext_str_out] ) ;
    
    % 检查输出文件是否存在
    if exist( file_out ) ~= 0 
        answer = input( sprintf( '文件 %s 已经存在，是否覆盖? ( Yes / [No] ):  ', file_out ), 's' ) ;
        if length( answer ) == 0
            answer = 'no' ;
        end
        
        answer = lower( answer ) ;
        if answer ~= 'y' | answer ~= 'yes' 
            disp( sprintf( '请使用另外的文件名，或备份已有的文件' ) ) ;
            disp( sprintf( '程序终止' ) );
            return ; 
        end
    end


    FemModel( file_in ) ;       % 定义有限元模型
    SolveModel ;                % 求解有限元模型
    SaveResults( file_out ) ;   % 保存计算结果
    DisplayResults(file_out) ;  % 显示结果
return ;

function FemModel( file_in )
%  定义有限元模型
%  输入参数：
%      file_in --- 输入文件名称
%  返回值：
%  说明：
    global gNode gElement gMaterial gBC1 gDF gNF
    
    % 打开文件
    fid = fopen( file_in, 'r' ) ;
    
    % 节点坐标
    node_number = fscanf( fid, '%d', 1 ) ;
    gNode = zeros( node_number, 3 ) ;
    for i = 1:node_number
        dummy = fscanf( fid, '%d', 1 ) ;
        gNode( i, : ) = fscanf( fid, '%f', [1 3] ) ;
    end
    
    % 单元定义
    element_number = fscanf( fid, '%d', 1 ) ;
    gElement = zeros( element_number, 5 ) ;
    for i = 1:element_number
        dummy = fscanf( fid, '%d', 1 ) ;
        gElement( i, : ) = fscanf( fid, '%d', [1 5] ) ;
    end
    
    % 材料性质
    material_number = fscanf( fid, '%d', 1 ) ;
    gMaterial = zeros( material_number, 2 ) ;
    for i = 1:material_number
        dummy = fscanf( fid, '%d', 1 ) ;
        gMaterial( i, : ) = fscanf( fid, '%f', [1 2] ) ;
    end
    
    % 第一类约束条件
    bc1_number = fscanf( fid, '%d', 1 ) ;
    gBC1 = zeros( bc1_number, 3 ) ;
    for i=1:bc1_number
        gBC1(i,1) = fscanf( fid, '%d', 1 ) ;  % 结点号
        gBC1(i,2) = fscanf( fid, '%d', 1 ) ;  % 自由度号
        gBC1(i,3) = fscanf( fid, '%f', 1 ) ;  % 约束值
    end
    
    % 分布压力（线性分布）
    df_number = fscanf( fid, '%d', 1 ) ;
    gDF = zeros( df_number, 7 ) ;
    for i=1:df_number
        gDF(i,1:3) = fscanf( fid, '%d', [1 3] ) ;  % 三个结点号
        gDF(i,4:6) = fscanf( fid, '%f', [1 3] ) ;  % 三个结点上的压力值 
    end
    
    % 集中力
    nf_number = fscanf( fid, '%d', 1 ) ;
    gNF = zeros( nf_number, 3 ) ;
    for i=1:nf_number
        gNF(i,1) = fscanf( fid, '%d', 1 ) ;  % 结点号
        gNF(i,2) = fscanf( fid, '%d', 1 ) ;  % 自由度号
        gNF(i,3) = fscanf( fid, '%f', 1 ) ;  % 集中力大小
    end
    
    % 关闭文件
    fclose( fid ) ;
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

    global gNode gElement gMaterial gBC1 gDF gNF gK gDelta gNodeStress gElementStress

    % step1. 定义整体刚度矩阵和节点力向量
    [node_number,dummy] = size( gNode ) ;
    gK = sparse( node_number * 3, node_number * 3 ) ;
    f = sparse( node_number * 3, 1 ) ;

    % step2. 计算单元刚度矩阵，并集成到整体刚度矩阵中
    [element_number,dummy] = size( gElement ) ;
    for ie=1:1:element_number
        disp( sprintf(  '计算节点刚度矩阵并集成，当前单元: %d', ie  ) ) ;
        k = StiffnessMatrix( ie ) ;
        AssembleStiffnessMatrix( ie, k ) ;
    end

    % step3. 计算分布压力的等效节点力，并集成到整体节点力向量中
    [df_number, dummy] = size( gDF ) ;
    for idf = 1:1:df_number
        enf = EquivalentDistPressure( gDF(idf,1:3), gDF(idf,4:6) ) ;
        i = gDF(idf, 1) ;
        j = gDF(idf, 2) ;
        m = gDF(idf, 3) ;
        f( (i-1)*3+1:(i-1)*3+3 ) = f( (i-1)*3+1:(i-1)*3+3 ) + enf( 1:3 ) ;
        f( (j-1)*3+1:(j-1)*3+3 ) = f( (j-1)*3+1:(j-1)*3+3 ) + enf( 4:6 ) ;
        f( (m-1)*3+1:(m-1)*3+3 ) = f( (m-1)*3+1:(m-1)*3+3 ) + enf( 7:9 ) ;
    end
  
    % step4. 把集中力集成到整体节点力向量中
    [nf_number, dummy] = size( gNF ); 
    for i = 1:1:nf_number
        in = gNF( i, 1 ) ;
        id = gNF( i, 2 ) ;
        f( (in-1)*3 + id ) = f( (in-1)*3 + id ) + gNF( i, 3 ) ;
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
    
    % step 6. 计算单元应力
    gElementStress = zeros( element_number, 6 ) ;
    for ie=1:element_number 
        disp( sprintf(  '计算单元应力，当前单元: %d', ie  ) ) ;
        es = ElementStress( ie ) ;
        gElementStress( ie, : ) = es ;
    end
    
    % step 7. 计算节点应力(采用绕节点加权平均)
    gNodeStress = zeros( node_number, 6 ) ;       
    for i=1:node_number                         
        disp( sprintf(  '计算节点应力，当前结点: %d', i  ) ) ;
        S = zeros( 1, 6 ) ;                         
        V = 0 ;
        for ie=1:1:element_number
            for k=1:1:4
                if i == gElement( ie, k ) 
                    vol = ElementVolume( ie ) ;
                    S = S + gElementStress(ie,:) * vol ;
                    V = V + vol ;
                    break ;
                end
            end
        end
        gNodeStress(i,:) = S / V ;
    end
return

function k = StiffnessMatrix( ie )
%  计算单元刚度矩阵
%  输入参数:
%     ie ----  单元号
%  返回值:
%     k  ----  单元刚度矩阵

    global gNode gElement gMaterial 
    
    % 读取结点坐标
    x = gNode( gElement( ie, : ), 1 ) ;
    y = gNode( gElement( ie, : ), 2 ) ;
    z = gNode( gElement( ie, : ), 3 ) ;
    
    % 计算 6V
    V6 = det( [1 x(1) y(1) z(1)
               1 x(2) y(2) z(2)
               1 x(3) y(3) z(3)
               1 x(4) y(4) z(4) ] ) ;
           
    if V6 < 0
        disp( sprintf( '警告：单元 %d 的结点排列顺序有问题', ie ) ) ;
        pause ;
    end
    % 计算应变矩阵
    B = zeros( 6, 12 ) ;
    for i=1:4
        j = mod(i,  4) + 1 ;
        m = mod(i+1,4) + 1 ;
        p = mod(i+2,4) + 1 ;
        bi = - det( [ 1 y(j) z(j)
                      1 y(m) z(m)
                      1 y(p) z(p) ] ) ;
        ci =   det( [ 1 x(j) z(j)
                      1 x(m) z(m)
                      1 x(p) z(p) ] ) ;
        di = - det( [ 1 x(j) y(j)
                      1 x(m) y(m)
                      1 x(p) y(p) ] ) ;
        B(:,(i-1)*3+1:(i-1)*3+3 ) = (-1)^(i+1)*[bi   0   0
                                                 0  ci   0
                                                 0   0  di
                                                 0  di  ci
                                                di   0  bi
                                                ci  bi   0] ;
    end
    B = B / V6 ;
    
    % 计算弹性矩阵
    E = gMaterial( gElement( ie, 5 ), 1 ) ;
    mu = gMaterial( gElement( ie, 5 ), 2 ) ;
    A1 = mu/(1-mu) ;
    A2 = (1-2*mu)/2/(1-mu) ;
    A3 = E*(1-mu)/(1+mu)/(1-2*mu) ;
    D = [  1  A1  A1   0   0   0
          A1   1  A1   0   0   0
          A1  A1   1   0   0   0
           0   0   0  A2   0   0 
           0   0   0  0   A2   0
           0   0   0  0    0  A2] ;
    D = D * A3 ;
    
    % 计算单元刚度矩阵
    V = abs(V6)/6 ;
    k = transpose( B ) * D * B * V ;
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
            for p=1:1:3
                for q=1:1:3
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

function enf = EquivalentDistPressure( node, pressure )
%   计算线性分布压力的等效节点力
%   输入参数:
%      node ----------  结点号
%      pressure ------  跟结点号对应的压力值
%   返回值:
%      enf -----------  等效节点力向量
    global gNode
    
    enf = zeros( 9, 1 ) ;
    % 计算作用压力的三角形的面积
    x = gNode( node, 1 ) ;
    y = gNode( node, 2 ) ;
    z = gNode( node, 3 ) ;
    b = - det( [1  y(1) z(1)
                1  y(2) z(2)
                1  y(3) z(3)] ) ;
    c =   det( [1  x(1) z(1)
                1  x(2) z(2)
                1  x(3) z(3)] ) ;
    d = - det( [1  x(1) y(1)
                1  x(2) y(2)
                1  x(3) y(3)] ) ;
    A2 = sqrt( b^2 + c^2 + d^2 ) ;
    A = A2 / 2 ;

    % 计算三个结点上的等效节点力
    f1 = ( pressure(1)   + pressure(2)/2 + pressure(3)/2 ) / 6 * A ;
    f2 = ( pressure(1)/2 + pressure(2)   + pressure(3)/2 ) / 6 * A ;
    f3 = ( pressure(1)/2 + pressure(2)/2 + pressure(3)   ) / 6 * A ;
    
    % 计算三个方向上的分量
    nx = -b/A2 ;
    ny = -c/A2 ;
    nz = -d/A2 ;
    enf(1:3) = f1*[nx;ny;nz] ;
    enf(4:6) = f2*[nx;ny;nz] ;
    enf(7:9) = f3*[nx;ny;nz] ;
return

function es = ElementStress( ie )
%  计算单元的应力分量
%  输入参数
%      ie  ----- 单元号
%  返回值
%      es ----- 单元应力分量列阵（1×6）： [sx, sy, sz, tyz, txz, txy]
    global gElement gDelta gMaterial
    
    de = zeros( 12, 1 ) ;   % 单元节点位移列阵
    E  = gMaterial( gElement(ie, 5), 1 ) ;
    mu = gMaterial( gElement(ie, 5), 2 ) ;
    A1 = mu/(1-mu) ;
    A2 = (1-2*mu)/2/(1-mu) ;
    A3 = E*(1-mu)/(1+mu)/(1-2*mu) ;
    D = [  1  A1  A1   0   0   0
          A1   1  A1   0   0   0
          A1  A1   1   0   0   0
           0   0   0  A2   0   0 
           0   0   0   0  A2   0
           0   0   0   0   0  A2] ;
    D = D * A3 ;
    B = MatrixB( ie ) ;
    for j=1:1:4
        de( 3*j-2 ) = gDelta( 3*gElement( ie, j )-2 ) ;
        de( 3*j-1 ) = gDelta( 3*gElement( ie, j )-1 ) ;
        de( 3*j   ) = gDelta( 3*gElement( ie, j )   ) ;
    end
    es = D * B * de ;
    es = transpose( es ) ;
return

function vol = ElementVolume( ie )
%  计算单元的体积
%  输入参数
%      ie  ----- 单元号
%  返回值
%      vol ----- 单元体积

    global gNode gElement
    
    % 读取结点坐标
    x = gNode( gElement( ie, : ), 1 ) ;
    y = gNode( gElement( ie, : ), 2 ) ;
    z = gNode( gElement( ie, : ), 3 ) ;
    
    % 计算 6V
    V6 = det( [1 x(1) y(1) z(1)
               1 x(2) y(2) z(2)
               1 x(3) y(3) z(3)
               1 x(4) y(4) z(4) ] ) ;
    
    vol = abs( V6/6 ) ;    
return

function B = MatrixB(ie)
%  计算单元的应变矩阵B
%  输入参数:
%     ie ----  单元号
%  返回值:
%     B  ----  单元应变矩阵

    global gNode gElement 
    
    % 读取结点坐标
    x = gNode( gElement( ie, : ), 1 ) ;
    y = gNode( gElement( ie, : ), 2 ) ;
    z = gNode( gElement( ie, : ), 3 ) ;
    
    % 计算 6V
    V6 = det( [1 x(1) y(1) z(1)
               1 x(2) y(2) z(2)
               1 x(3) y(3) z(3)
               1 x(4) y(4) z(4) ] ) ;
           
    % 计算应变矩阵
    B = zeros( 6, 12 ) ;
    for i=1:4
        j = mod(i,  4) + 1 ;
        m = mod(i+1,4) + 1 ;
        p = mod(i+2,4) + 1 ;
        bi = - det( [ 1 y(j) z(j)
                      1 y(m) z(m)
                      1 y(p) z(p) ] ) ;
        ci =   det( [ 1 x(j) z(j)
                      1 x(m) z(m)
                      1 x(p) z(p) ] ) ;
        di = - det( [ 1 x(j) y(j)
                      1 x(m) y(m)
                      1 x(p) y(p) ] ) ;
        B(:,(i-1)*3+1:(i-1)*3+3 ) = (-1)^(i+1)*[bi   0   0
                                                 0  ci   0
                                                 0   0  di
                                                 0  di  ci
                                                di   0  bi
                                                ci  bi   0] ;
    end
    B = B / V6 ;
return

function SaveResults( file_out ) 
%  保存计算结果
%  输入参数：
%     file_out ---- 存盘文件名
%  返回值：
%     无

    global gNode gElement gMaterial gBC1 gDF gNF gK gDelta gNodeStress gElementStress
    save( file_out, 'gNode', 'gElement', 'gMaterial', 'gBC1', ...
          'gDF', 'gNF', 'gDelta', 'gNodeStress', 'gElementStress' ) ;
return

function DisplayResults(file_out)
%  显示计算结果
%  输入参数：
%     无
%  返回值：
%     无
    y = [3,3.5,4,4.5] ;
    color = 'bgrc' ;
    [x,sy,syi] = exam5_1_post( file_out,y,0:0.25:4.5 ) ;
    figure ;
    hold ;
    for i=1:length(y)
        plot( x, sy(i,:),color(i) ) ;
    end
    legend( 'y=3.0','y=3.5','y=4.0','y=4.5' ) ;
    title( '离集中力不同距离处，正应力沿横截面的分布' ) ;
    xlabel( 'x' ) ;
    ylabel( '正应力' ) ;
    hold off ;
return