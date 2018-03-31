function exam5_2
% 本程序为第五章的第二个算例，采用三角形轴对称单元计算圆形垂直荷载作用下地基的变形

% 定义全局变量
%      gNode ------------- 节点坐标
%      gElement ---------- 单元定义
%      gMaterial --------- 材料性质
%      gBC1 -------------- 第一类约束条件
%      gK ---------------- 整体刚度矩阵
%      gDelta ------------ 整体节点坐标
    global gNode gElement gMaterial gBC1 gK gDelta

    file_in = 'exam5_2.dat' ;
    
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
    DisplayResults ;               % 显示计算结果    
return ;

function FemModel(filename)
%  定义有限元模型
%  输入参数：
%      filename --- 有限元模型文件
%  返回值：
%      无
%  说明：
%      该函数定义有限元模型数据：
%        gNode ------- 节点定义
%        gElement ---- 单元定义
%        gMaterial --- 材料定义，包括弹性模量，梁的截面积和梁的抗弯惯性矩
%        gBC1--------- 第一类约束条件
%        gDFStep ----- 分布力荷载步
%        gDF --------- 分布力

    global gNode gElement gMaterial gBC1 gDF gDFStep
    
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
    gMaterial = zeros( material_number, 2 ) ;
    for i=1:material_number
        dummy = fscanf( fid, '%d', 1 ) ;
        gMaterial( i, : ) = fscanf( fid, '%f', [1,2] ) ;
    end
    
    % 读取边界条件
    bc1_number = fscanf( fid, '%d', 1 ) ;
    gBC1 = zeros( bc1_number, 3 ) ;
    for i=1:bc1_number
        gBC1( i, 1 ) = fscanf( fid, '%d', 1 ) ;
        gBC1( i, 2 ) = fscanf( fid, '%d', 1 ) ;
        gBC1( i, 3 ) = fscanf( fid, '%f', 1 ) ;
    end
    
    % 读取分布压力步数
    dfstep_number = fscanf( fid, '%d', 1 ) ;
    gDFStep = fscanf( fid, '%d', [dfstep_number,1] ) ;
    
    % 读取分布压力
    gDF = [] ;
    for idfstep = 1:dfstep_number
        df_number = fscanf( fid, '%d', 1 ) ;
        DF = zeros( df_number, 4 ) ;
        for i=1:df_number
            DF( i, 1:2 ) = fscanf( fid, '%d', [1, 2] ) ;   % 结点号
            DF( i, 3:4 ) = fscanf( fid, '%f', [1, 2] ) ;   % 对应于结点号的压力值
        end
        gDF = [ gDF; DF ] ;
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

    global gNode gElement gMaterial gBC1 gDF gDFStep gK gDelta

    % step1. 定义整体刚度矩阵和节点力向量
    [node_number,dummy] = size( gNode ) ;
    gK = sparse( node_number * 2, node_number * 2 ) ;
    f = sparse( node_number * 2, length(gDFStep) ) ;

    % step2. 计算单元刚度矩阵，并集成到整体刚度矩阵中
    [element_number,dummy] = size( gElement ) ;
    for ie=1:1:element_number
        disp( sprintf(  '计算刚度矩阵，当前单元: %d', ie  ) ) ;
        k = StiffnessMatrix( ie ) ;
        AssembleStiffnessMatrix( ie, k ) ;
    end

    % step3. 计算分布压力产生的等效节点力
    for idfstep = 1:length(gDFStep)
        df_number = gDFStep( idfstep ) ;
        for i=1+sum(gDFStep(1:idfstep-1)):df_number+sum(gDFStep(1:idfstep-1))
            enf = EquivalentDistPressure( gDF(i,1), gDF(i,2), gDF(i,3), gDF(i,4) ) ;
            m = gDF(i, 1) ;
            n = gDF(i, 2) ;
            f(m*2-1:m*2,idfstep) = f(m*2-1:m*2,idfstep)+enf( 1:2 ) ;
            f(n*2-1:n*2,idfstep) = f(n*2-1:n*2,idfstep)+enf( 3:4 ) ;
        end
    end
    
    % step4. 处理约束条件，修改刚度矩阵和节点力向量。采用划行划列法
    [bc_number,dummy] = size( gBC1 ) ;
    for ibc=1:1:bc_number
        n = gBC1(ibc, 1 ) ;
        d = gBC1(ibc, 2 ) ;
        m = (n-1)*2 + d ;
        for idfstep = 1:length(gDFStep)
            f(:,idfstep) = f(:,idfstep) - gBC1(ibc,3) * gK(:,m) ;
            f(m,idfstep) = gBC1(ibc,3) ;
        end
        gK(:,m) = zeros( node_number*2, 1 ) ;
        gK(m,:) = zeros( 1, node_number*2 ) ;
        gK(m,m) = 1.0 ;
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
    E  = gMaterial( gElement(ie, 4), 1 ) ;
    mu = gMaterial( gElement(ie, 4), 2 ) ;
    ri = gNode( gElement( ie, 1 ), 1 ) ;
    zi = gNode( gElement( ie, 1 ), 2 ) ;
    rj = gNode( gElement( ie, 2 ), 1 ) ;
    zj = gNode( gElement( ie, 2 ), 2 ) ;
    rm = gNode( gElement( ie, 3 ), 1 ) ;
    zm = gNode( gElement( ie, 3 ), 2 ) ;
    ai = rj*zm - rm*zj ;
    aj = rm*zi - ri*zm ;
    am = ri*zj - rj*zi ;
    bi = zj - zm ;
    bj = zm - zi ;
    bm = zi - zj ;
    ci = -(rj-rm) ;
    cj = -(rm-ri) ;
    cm = -(ri-rj) ;
    area = (ai+aj+am)/2 ;
    r = (ri+rj+rm)/3; 
    z = (zi+zj+zm)/3 ;
    gi = ai/r + bi + ci*z/r ;
    gj = aj/r + bj + cj*z/r ;
    gm = am/r + bm + cm*z/r ;
    B = [bi  0 bj  0 bm  0
         gi  0 gj  0 gm  0
          0 ci  0 cj  0 cm
         ci bi cj bj cm bm] ;
    B = B/2/area ;
    A1 = mu/(1-mu) ;
    A2 = (1-2*mu)/2/(1-mu) ;
    A3 = E*(1-mu)/(1+mu)/(1-2*mu) ;
    D = A3*[ 1  A1  A1    0
            A1   1  A1    0
            A1  A1   1    0  
             0   0   0   A2] ;
    k = 2*pi*r*abs(area)*transpose(B)*D*B ;    
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

function enf = EquivalentDistPressure( m, n, pm, pn )
% 计算线性分布压力的等效节点力
% 输入参数:
%      i,j ----------  结点号
%      pi,pj --------  跟结点号对应的压力值
% 返回值:
%      enf -----------  等效节点力向量
    global gNode
    enf = zeros( 4, 1 ) ;
    
    % 获取结点坐标
    rm = gNode( m, 1 ) ;
    rn = gNode( n, 1 ) ;
    zm = gNode( m, 2 ) ;
    zn = gNode( n, 2 ) ;

    % 计算pi_, pj_
    pm_ = pi/6*( pm*(3*rm+rn) + pn*(rm+rn) ) ;
    pn_ = pi/6*( pm*(rm+rn) + pn*(rm+3*rn) ); 
    
    % 计算三个结点上的等效节点力
    enf(1:2) = pm_ * [zm-zn;rn-rm] ;
    enf(3:4) = pn_ * [zm-zn;rn-rm] ;
return

function SaveResults( file_out )
%  保存计算结果
%  输入参数：
%     file_out  --- 保存结果文件
%  返回值：
%     无
    global gNode gElement gMaterial gBC1 gDelta
    save( file_out, 'gNode', 'gElement', 'gMaterial', 'gBC1', 'gDelta' ) ;
return

function DisplayResults
%  显示计算结果
%  输入参数：
%     无
%  返回值：
%     无
    global gNode gElement gMaterial gBC1 gDelta
    
    %  解析解
    a = [ 1, 0.875, 0.750, 0.625, 0.5 ] ;
    width = 10 ;
    height = 15 ;
    q = 1e5 ;
    E = gMaterial( 1, 1 ) ;
    mu = gMaterial( 1, 2 ) ;
    w = -2*q*a*(1-mu^2)/E ;
    
    % 有限元解
    node = 47 ;   % 原点处的结点号
    wfem = full(gDelta( node*2, : )) ;

    % 显示结果
    fprintf( '\n\n\n               不同载荷作用半径下载荷圆心处的竖向位移\n' ) ;
    fprintf( '=======================================================================\n' ) ;
    fprintf( '载荷半径  计算区域(m)      有限元        解析解        相对误差\n' ) ;
    fprintf( '  (m)    (深度×宽度)       (mm)          (mm)             (%%) \n' ) ;
    fprintf( '-----------------------------------------------------------------------\n' ) ;
    for i=1:length(a)
        fprintf( '%6.3f     %2.0f×%2.0f       %8.3f       %8.3f        %8.3f\n', ...
                  a(i), height, width, wfem(i)*1000, w(i)*1000, (wfem(i)-w(i))/w(i)*100 ) ;
    end
    fprintf( '=======================================================================\n' ) ;
return

