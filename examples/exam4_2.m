function exam4_2( file_in )
% 本程序为第四章的第二个算例，采用三角形单元计算隧道在自重作用下的变形和应力
%      exam4_2( filename )
%  输入参数： 
%      file_in  ---------- 有限元模型文件

% 定义全局变量
%      gNode ------------- 节点坐标
%      gElement ---------- 单元定义
%      gMaterial --------- 材料性质
%      gBC1 -------------- 第一类约束条件
%      gK ---------------- 整体刚度矩阵
%      gDelta ------------ 整体节点坐标
%      gNodeStress ------- 节点应力
%      gElementStress ---- 单元应力
    global gNode gElement gMaterial gBC1 gK gDelta gNodeStress gElementStress

    if nargin < 1
        file_in = 'exam4_2.dat' ;
    end
    
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

    % 建立有限元模型并求解，保存结果
    FemModel( file_in ) ;          % 定义有限元模型
    SolveModel ;                   % 求解有限元模型
    SaveResults( file_out ) ;      % 保存计算结果
    
    % 计算结束
    disp( sprintf( '计算正常结束，结果保存在文件 %s 中', file_out ) ) ;
    disp( sprintf( '可以使用后处理程序 exam4_2_post.m 显示计算结果' ) ) ;
return ;

function FemModel(filename)
%  定义有限元模型
%  输入参数：
%      filename --- 有限元模型文件
%  返回值：
%      无
%  说明：
%      该函数定义平面问题的有限元模型数据：
%        gNode ------- 节点定义
%        gElement ---- 单元定义
%        gMaterial --- 材料定义，包括弹性模量，梁的截面积和梁的抗弯惯性矩
%        gBC1 -------- 约束条件

    global gNode gElement gMaterial gBC1
    
    % 打开文件
    fid = fopen( filename, 'r' ) ;
    
    % 读取节点坐标
    node_number = fscanf( fid, '%d', 1 ) ;
    gNode = zeros( node_number, 2 ) ;
    for i=1:node_number
        dummy = fscanf( fid, '%d', 1 ) ;
        gNode( i, : ) = fscanf( fid, '%f', [1, 2] ) ;
    end
    
    % 读取单元定义
    element_number = fscanf( fid, '%d', 1 ) ;
    gElement = zeros( element_number, 4 ) ;
    for i=1:element_number
        dummy = fscanf( fid, '%d', 1 ) ;
        gElement( i, : ) = fscanf( fid, '%d', [1, 4] ) ;
    end
    
    % 读取材料信息
    material_number = fscanf( fid, '%d', 1 ) ;
    gMaterial = zeros( material_number, 4 ) ;
    for i=1:material_number
        dummy = fscanf( fid, '%d', 1 ) ;
        gMaterial( i, : ) = fscanf( fid, '%f', [1,4] ) ;
    end
    
    % 读取边界条件
    bc1_number = fscanf( fid, '%d', 1 ) ;
    gBC1 = zeros( bc1_number, 3 ) ;
    for i=1:bc1_number
        gBC1( i, 1 ) = fscanf( fid, '%d', 1 ) ;
        gBC1( i, 2 ) = fscanf( fid, '%d', 1 ) ;
        gBC1( i, 3 ) = fscanf( fid, '%f', 1 ) ;
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
%        5. 计算单元应力和节点应力

    global gNode gElement gMaterial gBC1 gK gDelta gNodeStress gElementStress

    % step1. 定义整体刚度矩阵和节点力向量
    [node_number,dummy] = size( gNode ) ;
    gK = sparse( node_number * 2, node_number * 2 ) ;
    f = sparse( node_number * 2, 1 ) ;

    % step2. 计算单元刚度矩阵，并集成到整体刚度矩阵中
    [element_number,dummy] = size( gElement ) ;
    for ie=1:1:element_number
        disp( sprintf(  '计算刚度矩阵，当前单元: %d', ie  ) ) ;
        k = StiffnessMatrix( ie ) ;
        AssembleStiffnessMatrix( ie, k ) ;
    end

    % step3. 计算自重产生的等效节点力
    for ie=1:1:element_number
        disp( sprintf(  '计算自重的等效节点力，当前单元: %d', ie  ) ) ;
        egf = EquivalentGravityForce( ie ) ;
        i = gElement( ie, 1 ) ;
        j = gElement( ie, 2 ) ;
        m = gElement( ie, 3 ) ;
        f( (i-1)*2+1 : (i-1)*2+2 ) = f( (i-1)*2+1 : (i-1)*2+2 ) + egf( 1:2 ) ;
        f( (j-1)*2+1 : (j-1)*2+2 ) = f( (j-1)*2+1 : (j-1)*2+2 ) + egf( 3:4 ) ;
        f( (m-1)*2+1 : (m-1)*2+2 ) = f( (m-1)*2+1 : (m-1)*2+2 ) + egf( 5:6 ) ;
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
        S = zeros( 1, 3 ) ;                         
        A = 0 ;
        for ie=1:1:element_number
            for k=1:1:3
                if i == gElement( ie, k ) 
                    area= ElementArea( ie ) ;
                    S = S + gElementStress(ie,1:3 ) * area ;
                    A = A + area ;
                    break ;
                end
            end
        end
        gNodeStress(i,1:3) = S / A ;
        gNodeStress(i,6) = 0.5*sqrt( (gNodeStress(i,1)-gNodeStress(i,2))^2 + 4*gNodeStress(i,3)^2 ) ;
        gNodeStress(i,4) = 0.5*(gNodeStress(i,1)+gNodeStress(i,2)) + gNodeStress(i,6) ;
        gNodeStress(i,5) = 0.5*(gNodeStress(i,1)+gNodeStress(i,2)) - gNodeStress(i,6) ;
    end
return

function k = StiffnessMatrix( ie )
%  计算单元刚度矩阵
%  输入参数:
%     ie ----  单元号
%  返回值:
%     k  ----  单元刚度矩阵

    global gNode gElement gMaterial 
    k = zeros( 6, 6 ) ;
    E  = gMaterial( gElement(ie, 4), 1 ) ;
    mu = gMaterial( gElement(ie, 4), 2 ) ;
    h  = gMaterial( gElement(ie, 4), 3 ) ;
    xi = gNode( gElement( ie, 1 ), 1 ) ;
    yi = gNode( gElement( ie, 1 ), 2 ) ;
    xj = gNode( gElement( ie, 2 ), 1 ) ;
    yj = gNode( gElement( ie, 2 ), 2 ) ;
    xm = gNode( gElement( ie, 3 ), 1 ) ;
    ym = gNode( gElement( ie, 3 ), 2 ) ;
    ai = xj*ym - xm*yj ;
    aj = xm*yi - xi*ym ;
    am = xi*yj - xj*yi ;
    bi = yj - ym ;
    bj = ym - yi ;
    bm = yi - yj ;
    ci = -(xj-xm) ;
    cj = -(xm-xi) ;
    cm = -(xi-xj) ;
    area = (ai+aj+am)/2 ;
    B = [bi  0 bj  0 bm  0
          0 ci  0 cj  0 cm
         ci bi cj bj cm bm] ;
    B = B/2/area ;
    D = [ 1-mu    mu      0
           mu    1-mu     0
            0      0   (1-2*mu)/2] ;
    D = D*E/(1-2*mu)/(1+mu) ;
    k = transpose(B)*D*B*h*abs(area) ;    
return

function B = MatrixB( ie )
%  计算单元的应变矩阵B
%  输入参数:
%     ie ----  单元号
%  返回值:
%     B  ----  单元应变矩阵
    global gNode gElement
    xi = gNode( gElement( ie, 1 ), 1 ) ;
    yi = gNode( gElement( ie, 1 ), 2 ) ;
    xj = gNode( gElement( ie, 2 ), 1 ) ;
    yj = gNode( gElement( ie, 2 ), 2 ) ;
    xm = gNode( gElement( ie, 3 ), 1 ) ;
    ym = gNode( gElement( ie, 3 ), 2 ) ;
    ai = xj*ym - xm*yj ;
    aj = xm*yi - xi*ym ;
    am = xi*yj - xj*yi ;
    bi = yj - ym ;
    bj = ym - yi ;
    bm = yi - yj ;
    ci = -(xj-xm) ;
    cj = -(xm-xi) ;
    cm = -(xi-xj) ;
    area = (ai+aj+am)/2 ;
    B = [bi  0 bj  0 bm  0
          0 ci  0 cj  0 cm
         ci bi cj bj cm bm] ;
    B = B/2/area ;
return

function area = ElementArea( ie )
%  计算单元面积
%  输入参数:
%     ie ----  单元号
%  返回值:
%     area  ----  单元面积
    global gNode gElement gMaterial 

    xi = gNode( gElement( ie, 1 ), 1 ) ;
    yi = gNode( gElement( ie, 1 ), 2 ) ;
    xj = gNode( gElement( ie, 2 ), 1 ) ;
    yj = gNode( gElement( ie, 2 ), 2 ) ;
    xm = gNode( gElement( ie, 3 ), 1 ) ;
    ym = gNode( gElement( ie, 3 ), 2 ) ;
    ai = xj*ym - xm*yj ;
    aj = xm*yi - xi*ym ;
    am = xi*yj - xj*yi ;
    area = abs((ai+aj+am)/2) ;
return

function AssembleStiffnessMatrix( ie, k )
%  把单元刚度矩阵集成到整体刚度矩阵
%  输入参数:
%      ie  --- 单元号
%      k   --- 单元刚度矩阵
%  返回值:
%      无
    global gElement gK
    for i=1:1:3
        for j=1:1:3
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

function egf = EquivalentGravityForce( ie )
%  计算单元自重的等效节点力
%  输入参数
%      ie  ----- 单元号
%  返回值
%      egf ----- 自重的等效节点力 
    global gElement gMaterial
    
    area = ElementArea( ie ) ;
    h    = gMaterial( gElement( ie, 4 ), 3 ) ;
    ro   = gMaterial( gElement( ie, 4 ), 4 ) ;
    w    = area * h * ro ;
    egf  = -w/3 * [0; 1; 0; 1; 0; 1] ;
return

function es = ElementStress( ie )
%  计算单元的应力分量
%  输入参数
%      ie  ----- 单元号
%  返回值
%      es ----- 单元应力分量列阵（1×6）： [sx, sy, txy, s1, s2, tmax]
    global gElement gDelta  gMaterial
    
    es = zeros( 1, 6 ) ;   % 单元的应力分量
    de = zeros( 6, 1 ) ;   % 单元节点位移列阵
    E  = gMaterial( gElement(ie, 4), 1 ) ;
    mu = gMaterial( gElement(ie, 4), 2 ) ;
    D = [ 1-mu    mu      0
           mu    1-mu     0
            0      0   (1-2*mu)/2] ;
    D = D*E/(1-2*mu)/(1+mu) ;
    B = MatrixB( ie ) ;
    for j=1:1:3
        de( 2*j-1 ) = gDelta( 2*gElement( ie, j )-1 ) ;
        de( 2*j   ) = gDelta( 2*gElement( ie, j )   ) ;
    end
    es(1:3) = D * B * de ;
    es(6) = 0.5*sqrt((es(1)-es(2))^2 + 4*es(3)^2 ) ;
    es(4) = 0.5*(es(1)+es(2)) + es(6) ;
    es(5) = 0.5*(es(1)+es(2)) - es(6) ;
return

function SaveResults( file_out )
%  显示计算结果
%  输入参数：
%     file_out  --- 保存结果文件
%  返回值：
%     无

    global gNode gElement gMaterial gBC1 gDelta gNodeStress gElementStress
    save( file_out, 'gNode', 'gElement', 'gMaterial', 'gBC1', 'gDelta', 'gNodeStress', 'gElementStress' ) ;
return