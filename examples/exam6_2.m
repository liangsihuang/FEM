function exam6_2
% 本程序为第六章的第二个算例，采用20结点六面体等参元计算端部受横向力的悬臂梁
%  输入参数： 
%       无

    file_in = 'exam6_2.dat' ;
    
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
    DisplayResults ;            % 把计算结果显示出来 
return

function FemModel( file_in )
%  定义有限元模型
%  输入参数：
%      file_in ---- 有限模型数据文件
%  返回值：
%      无

%  全局变量及名称
%      gNode ------ 节点坐标
%      gElement --- 单元定义
%      gMaterial -- 材料性质
%      gBC1 ------- 第一类约束条件
%      gNF -------- 集中力
    global gNode gElement gMaterial gBC1 gNF
    
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
    gElement = zeros( element_number, 21 ) ;
    for i = 1:element_number
        dummy = fscanf( fid, '%d', 1 ) ;
        gElement( i, : ) = fscanf( fid, '%d', [1 21] ) ;
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

%  全局变量及名称
%      gNode ----------- 节点坐标
%      gElement -------- 单元定义
%      gMaterial ------- 材料性质
%      gBC1 ------------ 第一类约束条件
%      gNF ------------- 集中力
%      gK -------------- 整体刚度矩阵
%      gDelta ---------- 结点位移向量
%      gElementStress -- 每个单元的结点应力
    global gNode gElement gMaterial gBC1 gNF gK gDelta gElementStress

    % step1. 定义整体刚度矩阵和节点力向量
    [node_number,dummy] = size( gNode ) ;
    %gK = sparse( node_number * 3, node_number * 3 ) ;   % 由于稀疏矩阵的存储速度很慢，因此改成普通矩阵运算
    %f = sparse( node_number * 3, 1 ) ;                  % 如果计算机内存不够，就恢复这两句，并注释掉下面两句。
    gK = zeros( node_number*3,node_number*3) ;
     f = zeros( node_number*3,1 ) ;
     
    % step2. 计算单元刚度矩阵，并集成到整体刚度矩阵中
    [element_number,dummy] = size( gElement ) ;
    for ie=1:1:element_number
        disp( sprintf(  '计算单元刚度矩阵并集成，当前单元: %d', ie  ) ) ;
        k = StiffnessMatrix( ie ) ;
        AssembleStiffnessMatrix( ie, k ) ;
    end
    
    spy(gK) ;
    nz=nnz(gK) ;
    fprintf( 'non-zeros elements/total elmements = %g%%',nz/node_number/node_number/9*100) ;
    pause ;
    
    % step3. 把集中力集成到整体节点力向量中
    [nf_number, dummy] = size( gNF ) ;
    for idf = 1:1:nf_number
        node  = gNF(idf, 1) ;
        dim   = gNF(idf, 2) ;
        force = gNF(idf, 3) ;
        f( (node-1)*3+dim ) = f( (node-1)*3+dim ) + force ;
    end
    
    % step4. 处理第一类约束条件，修改刚度矩阵和节点力向量。采用乘大数法
    [bc1_number,dummy] = size( gBC1 ) ;
    for ibc=1:1:bc1_number
        n = gBC1(ibc, 1 ) ;
        d = gBC1(ibc, 2 ) ;
        m = (n-1)*3 + d ;
        f(m,:) = gBC1(ibc, 3)* gK(m,m) * 1e15 ;
        gK(m,m) = gK(m,m) * 1e15 ;
    end
    
    % step5. 求解方程组，得到节点位移向量
    gDelta = gK \ f ;
    
    % step6. 计算每个单元的结点应力
    gElementStress = zeros( element_number, 20, 6 ) ;
    delta = zeros( 60, 1 ) ;
    x = [-1, 1, 1, -1, -1, 1, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, -1, 1, 1, -1];
    e = [-1, -1, 1, 1, -1, -1, 1, 1, -1, 1, 1, -1, 0, 0, 0, 0, -1, -1, 1, 1];
    z = [-1, -1, -1, -1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, 1, -1, 0, 0, 0, 0];
    for ie = 1:element_number
        fprintf( '计算单元应力，当前单元号：%d\n', ie ) ;
        for n=1:20
            B = MatrixB( ie, x(n), e(n), z(n) ) ;
            D = MatrixD( ie ) ;
            delta(1:3:58) = gDelta( gElement(ie,1:20)*3-2) ;
            delta(2:3:59) = gDelta( gElement(ie,1:20)*3-1) ;
            delta(3:3:60) = gDelta( gElement(ie,1:20)*3) ;
            sigma = D*B*delta ;
            gElementStress( ie, n, : ) = sigma ;
        end
    end
return

function k = StiffnessMatrix( ie )   
%  计算平面应变等参数单元的刚度矩阵
%  输入参数：
%     ie -- 单元号
%  返回值：
%     k  -- 单元刚度矩阵
%  说明：
%     用高斯积分法求解平面等参数单元的刚度矩阵

    k = zeros( 60, 60 ) ;
    D = MatrixD( ie ) ;
    %3 x 3 x 3 高斯积分点和权系数
    %x = [-0.774596669241483,              0.0,0.774596669241483] ;   
    %w = [ 0.555555555555556,0.888888888888889,0.555555555555556] ;  
    %for i=1:1:length(x)
    %    for j=1:1:length(x)
    %        for m=1:1:length(x)
    %            B = MatrixB( ie, x(i), x(j), x(m) ) ;
    %            J = Jacobi( ie, x(i), x(j), x(m) ) ;
    %            k = k + w(i)*w(j)*w(m)*transpose(B)*D*B*det(J) ;   
    %        end
    %    end
    %end
    
    %14点 irons 积分, 详细参见以下文献 
    % Irons B. M., 1971 Quadrature rules for brick based finite elements, 
    %   International Journal for Numerical Methods in Engineering, 
    %   3(2):293-294.
    b0 = 320/361 ;
    c0 = 121/361 ;
    b  = sqrt( 19/30 ) ;
    c  = sqrt( 19/33 ) ;
    x1 = [-b,  0,  0
           b,  0,  0
           0, -b,  0
           0,  b,  0
           0,  0, -b
           0,  0,  b ] ;
    x2 = [-c, -c, -c
           c,  c,  c
           c, -c, -c
          -c,  c, -c
          -c, -c,  c
           c,  c, -c
          -c,  c,  c
           c, -c,  c ] ;
    for i=1:6
        B = MatrixB( ie, x1(i,1), x1(i,2), x1(i,3) ) ;
        J = Jacobi( ie, x1(i,1), x1(i,2), x1(i,3) ) ;
        k = k + b0*transpose(B)*D*B*det(J) ;
    end
    for i=1:8
        B = MatrixB( ie, x2(i,1), x2(i,2), x2(i,3) ) ;
        J = Jacobi( ie, x2(i,1), x2(i,2), x2(i,3) ) ;
        k = k + c0*transpose(B)*D*B*det(J) ;
    end
return

function D = MatrixD( ie )
%  计算单元的弹性矩阵D
%  输入参数：
%     ie --------- 单元号
%  返回值：
%     D  --------- 弹性矩阵D

%  全局变量及名称
%      gElement -------- 单元定义
%      gMaterial ------- 材料性质
    global gElement gMaterial 

    E  = gMaterial( gElement(ie, 21), 1 ) ;   % 弹性模量
    mu = gMaterial( gElement(ie, 21), 2 ) ;   % 泊松比 
    A1 = mu/(1-mu) ;                        
    A2 = (1-2*mu)/2/(1-mu) ;                
    A3 = E*(1-mu)/(1+mu)/(1-2*mu) ;         
    D = A3*[  1  A1  A1   0   0   0
             A1   1  A1   0   0   0
             A1  A1   1   0   0   0
              0   0   0  A2   0   0
              0   0   0   0  A2   0
              0   0   0   0   0  A2] ;
return

function B = MatrixB( ie, xi, eta, zeta )
%  计算单元的应变矩阵B
%  输入参数：
%     ie --------------- 单元号
%     xi,eta, zeta ----- 局部坐标  
%  返回值：
%     B  --------------- 在局部坐标处的应变矩阵B

    [N_x, N_y, N_z] = N_xyz( ie, xi, eta, zeta ); 
    B = zeros( 6, 60 ) ;
    for i=1:1:20
        B(:,(3*i-2):3*i) = [ N_x(i)       0       0
                                  0  N_y(i)       0
                                  0       0  N_z(i)
                                  0  N_z(i)  N_y(i)
                             N_z(i)      0   N_x(i)
                             N_y(i)  N_x(i)       0 ] ;
    end
return

function [N_x, N_y, N_z] = N_xyz( ie, xi, eta, zeta )
%  计算形函数对整体坐标的导数
%  输入参数：
%     ie -------------- 单元号
%     xi,eta,zeta ----- 局部坐标  
%  返回值：
%     N_x  ------------ 在局部坐标处的形函数对x的导数
%     N_y  ------------ 在局部坐标处的形函数对y的导数
%     N_z  ------------ 在局部坐标处的形函数对z的导数

    J = Jacobi( ie, xi, eta, zeta ) ;
    [N_xi,N_eta,N_zeta] = N_xietazeta( ie, xi, eta, zeta ) ;
    N = [N_xi;N_eta;N_zeta] ;
    A = inv(J)*N ;
    N_x = A(1,:) ;
    N_y = A(2,:) ;
    N_z = A(3,:) ;
return

function [N_xi, N_eta, N_zeta] = N_xietazeta( ie, xi, eta, zeta )
%  计算形函数对局部坐标的导数
%  输入参数：
%     ie -------------- 单元号
%     xi,eta,zeta ----- 局部坐标  
%  返回值：
%     N_xi   ---------- 在局部坐标处的形函数对xi的导数
%     N_eta  ---------- 在局部坐标处的形函数对eta的导数
%     N_zeta  --------- 在局部坐标处的形函数对zeta的导数

    x = [-1, 1, 1, -1, -1, 1, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, -1, 1, 1, -1];
    e = [-1, -1, 1, 1, -1, -1, 1, 1, -1, 1, 1, -1, 0, 0, 0, 0, -1, -1, 1, 1];
    z = [-1, -1, -1, -1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, 1, -1, 0, 0, 0, 0];
    N_xi   = zeros( 1, 20 ) ;
    N_eta  = zeros( 1, 20 ) ;
    N_zeta = zeros( 1, 20 ) ;

    N_xi(1:8)   = x(1:8).*(1+e(1:8)*eta).*(1+z(1:8)*zeta)/8 ;
    N_eta(1:8)  = (1+x(1:8)*xi).*e(1:8).*(1+z(1:8)*zeta)/8 ;
    N_zeta(1:8) = (1+x(1:8)*xi).*(1+e(1:8)*eta).*z(1:8)/8 ;

    N_xi(9:12)   = -2*xi*(1+e(9:12)*eta).*(1+z(9:12)*zeta)/4 ;
    N_eta(9:12)  = ( 1-xi^2)*e(9:12).*(1+z(9:12)*zeta)/4 ;
    N_zeta(9:12) = ( 1-xi^2)*(1+e(9:12)*eta).*z(9:12)/4 ;

    N_xi(13:16) = ( 1-eta^2)*x(13:16).*(1+z(13:16)*zeta)/4 ;
    N_eta(13:16) = -2*eta*(1+x(13:16)*xi).*(1+z(13:16)*zeta)/4 ;
    N_zeta(13:16) = ( 1-eta^2)*(1+x(13:16)*xi).*z(13:16)/4 ;

    N_xi(17:20) = ( 1-zeta^2)*x(17:20).*(1+e(17:20)*eta)/4 ;
    N_eta(17:20) = ( 1-zeta^2)*(1+x(17:20)*xi).*e(17:20)/4 ;
    N_zeta(17:20) = -2*zeta*(1+x(17:20)*xi).*(1+e(17:20)*eta)/4 ;

    N_xi(1) = N_xi(1) - 0.5*(N_xi(9) +N_xi(16)+N_xi(17)) ;
    N_xi(2) = N_xi(2) - 0.5*(N_xi(9) +N_xi(13)+N_xi(18)) ;
    N_xi(3) = N_xi(3) - 0.5*(N_xi(10)+N_xi(13)+N_xi(19)) ;
    N_xi(4) = N_xi(4) - 0.5*(N_xi(10)+N_xi(16)+N_xi(20)) ;
    N_xi(5) = N_xi(5) - 0.5*(N_xi(12)+N_xi(15)+N_xi(17)) ;
    N_xi(6) = N_xi(6) - 0.5*(N_xi(12)+N_xi(14)+N_xi(18)) ;
    N_xi(7) = N_xi(7) - 0.5*(N_xi(11)+N_xi(14)+N_xi(19)) ;
    N_xi(8) = N_xi(8) - 0.5*(N_xi(11)+N_xi(15)+N_xi(20)) ;

    N_eta(1) = N_eta(1) - 0.5*(N_eta(9) +N_eta(16)+N_eta(17)) ;
    N_eta(2) = N_eta(2) - 0.5*(N_eta(9) +N_eta(13)+N_eta(18)) ;
    N_eta(3) = N_eta(3) - 0.5*(N_eta(10)+N_eta(13)+N_eta(19)) ;
    N_eta(4) = N_eta(4) - 0.5*(N_eta(10)+N_eta(16)+N_eta(20)) ;
    N_eta(5) = N_eta(5) - 0.5*(N_eta(12)+N_eta(15)+N_eta(17)) ;
    N_eta(6) = N_eta(6) - 0.5*(N_eta(12)+N_eta(14)+N_eta(18)) ;
    N_eta(7) = N_eta(7) - 0.5*(N_eta(11)+N_eta(14)+N_eta(19)) ;
    N_eta(8) = N_eta(8) - 0.5*(N_eta(11)+N_eta(15)+N_eta(20)) ;

    N_zeta(1) = N_zeta(1) - 0.5*(N_zeta(9) +N_zeta(16)+N_zeta(17)) ;
    N_zeta(2) = N_zeta(2) - 0.5*(N_zeta(9) +N_zeta(13)+N_zeta(18)) ;
    N_zeta(3) = N_zeta(3) - 0.5*(N_zeta(10)+N_zeta(13)+N_zeta(19)) ;
    N_zeta(4) = N_zeta(4) - 0.5*(N_zeta(10)+N_zeta(16)+N_zeta(20)) ;
    N_zeta(5) = N_zeta(5) - 0.5*(N_zeta(12)+N_zeta(15)+N_zeta(17)) ;
    N_zeta(6) = N_zeta(6) - 0.5*(N_zeta(12)+N_zeta(14)+N_zeta(18)) ;
    N_zeta(7) = N_zeta(7) - 0.5*(N_zeta(11)+N_zeta(14)+N_zeta(19)) ;
    N_zeta(8) = N_zeta(8) - 0.5*(N_zeta(11)+N_zeta(15)+N_zeta(20)) ;
return

function J = Jacobi( ie, xi, eta, zeta )
%  计算Jacobi矩阵
%  输入参数：
%     ie -------------- 单元号
%     xi,eta,zeta ----- 局部坐标  
%  返回值：
%     J   ------------- 在局部坐标(xi,eta,zeta)处的Jacobi矩阵

%  全局变量及名称
%      gNode ----------- 节点坐标
%      gElement -------- 单元定义
    global gNode gElement

    x = gNode(gElement(ie,1:20),1) ;
    y = gNode(gElement(ie,1:20),2) ;
    z = gNode(gElement(ie,1:20),3) ;
    [N_xi, N_eta, N_zeta] = N_xietazeta( ie, xi, eta, zeta ) ;
    x_xi  = N_xi  * x ;
    x_eta = N_eta * x ;
    x_zeta = N_zeta * x ;

    y_xi  = N_xi  * y ;
    y_eta = N_eta * y ;
    y_zeta = N_zeta* y;

    z_xi = N_xi * z ;
    z_eta = N_eta * z ;
    z_zeta = N_zeta * z ;

    J = [  x_xi,    y_xi,    z_xi
          x_eta,   y_eta,   z_eta
         x_zeta,  y_zeta,  z_zeta ];
return

function AssembleStiffnessMatrix( ie, k )
%  把单元刚度矩阵集成到整体刚度矩阵
%  输入参数:
%      ie  --- 单元号
%      k   --- 单元刚度矩阵
%  返回值:
%      无

%  全局变量及名称
%      gElement -------- 单元定义
%      gK -------------- 整体刚度矩阵
    global gElement gK
    
%     for i=1:1:20
%         for j=1:1:20
%             for p=1:1:3
%                 for q=1:1:3
%                     m = (i-1)*3+p ;
%                     n = (j-1)*3+q ;
%                     M = (gElement(ie,i)-1)*3+p ;
%                     N = (gElement(ie,j)-1)*3+q ;
%                     gK(M,N) = gK(M,N) + k(m,n) ;
%                 end
%             end
%         end
%     end
    m = [(gElement(ie,1:20)-1)*3+1;
         (gElement(ie,1:20)-1)*3+2;
         (gElement(ie,1:20)-1)*3+3] ;
    gK(m(:),m(:)) = gK(m(:),m(:)) + k ;
return

function edf = EquivalentDistPressure( node, pressure )
%   计算分布压力的等效节点力
%   输入参数:
%      node ----------  结点号
%      pressure ------  跟结点号对应的压力值
%   返回值:
%      edf -----------  等效节点力向量

%  全局变量及名称
%      gNode ----------- 节点坐标
    global gNode
    
    x = gNode( node, 1 ) ;
    y = gNode( node, 2 ) ;
    z = gNode( node, 3 ) ;
    g = [-0.774596669241483,              0.0,0.774596669241483] ;   
    w = [ 0.555555555555556,0.888888888888889,0.555555555555556] ;  
    edf = zeros( 8*3, 1 ) ;
    for i=1:length(g)
        for j=1:length(g)
            xi = g(i) ;
            eta = g(j) ;
            N5 = ( eta - 1 ) * ( xi^2 - 1 ) / 2 ;
            N6 = ( xi + 1 ) * ( 1 - eta^2 ) / 2 ;
            N7 = ( eta + 1 ) * ( 1 - xi^2 ) / 2 ;
            N8 = ( xi - 1 ) * ( eta^2 - 1 ) / 2 ;
            N1 = ( 1 - xi ) * ( 1 - eta ) / 4 - 0.5 * ( N8 + N5 ) ;
            N2 = ( 1 + xi ) * ( 1 - eta ) / 4 - 0.5 * ( N5 + N6 ) ;
            N3 = ( 1 + xi ) * ( 1 + eta ) / 4 - 0.5 * ( N6 + N7 ) ;
            N4 = ( 1 - xi ) * ( 1 + eta ) / 4 - 0.5 * ( N7 + N8 ) ;
            N = [ N1 N2 N3 N4 N5 N6 N7 N8 ] ;
            N5_xi = (eta-1)*xi ;
            N6_xi = 1/2*(1-eta^2) ;
            N7_xi = -(eta+1)*xi ;
            N8_xi = 1/2*(eta^2-1) ;
            N1_xi = -1/4*(1-eta) - 1/2 * ( N5_xi + N8_xi) ;
            N2_xi = 1/4*(1-eta) - 1/2 * (N5_xi + N6_xi) ;
            N3_xi = 1/4*(eta+1) - 1/2 * (N6_xi + N7_xi) ;
            N4_xi = -1/4*(eta+1) - 1/2 * (N7_xi + N8_xi ) ;
            N_xi = [ N1_xi N2_xi N3_xi N4_xi N5_xi N6_xi N7_xi N8_xi ] ;
            N5_eta = 1/2*(xi^2-1) ;
            N6_eta = -eta*(1+xi) ;
            N7_eta = 1/2*(1-xi^2) ;
            N8_eta = eta*(xi-1) ;
            N1_eta = -1/4*(1-xi) - 1/2 * ( N5_eta + N8_eta) ;
            N2_eta = -1/4*(1+xi) - 1/2 * (N5_eta + N6_eta) ;
            N3_eta = 1/4*(1+xi) - 1/2 * (N6_eta + N7_eta) ;
            N4_eta = 1/4*(1-xi) - 1/2 * (N7_eta + N8_eta) ;
            N_eta = [ N1_eta N2_eta N3_eta N4_eta N5_eta N6_eta N7_eta N8_eta ] ;
            x_xi = N_xi * x ;
            y_xi = N_xi * y ;
            z_xi = N_xi * z ;
            x_eta = N_eta * x ;
            y_eta = N_eta * y ;
            z_eta = N_eta * z ;
            px = y_xi * z_eta - y_eta * z_xi ;
            py = z_xi * x_eta - z_eta * x_xi ;
            pz = x_xi * y_eta - x_eta * y_xi ;
            sig = N * transpose( pressure ) ;
            edf = edf + w(i)*w(j)*sig*[N1*eye(3);N2*eye(3);N3*eye(3); ...
              N4*eye(3);N5*eye(3);N6*eye(3);N7*eye(3);N8*eye(3)]*[px;py;pz] ;
        end
    end
return

function SaveResults( file_out ) 
%  保存计算结果
%  输入参数：
%     无
%  返回值：
%     无

%  全局变量及名称
%      gNode ----------- 节点坐标
%      gElement -------- 单元定义
%      gMaterial ------- 材料性质
%      gBC1 ------------ 第一类约束条件
%      gNF ------------- 集中力
%      gK -------------- 整体刚度矩阵
%      gDelta ---------- 结点位移向量
%      gElementStress -- 每个单元的结点应力
    global gNode gElement gMaterial gBC1 gDF gK gDelta gElementStress
    save( file_out, 'gNode', 'gElement', 'gMaterial', 'gBC1', ...
          'gDF', 'gDelta', 'gElementStress' ) ;
return

function DisplayResults
%  显示计算结果
%  输入参数：
%     无
%  返回值：
%     无
   % 读取正应力沿截面高度的变化
   global gElement gNode gMaterial gElementStress
   
    [element_number,dummy] = size( gElement ) ;
    node = [ 112  426  445  463  482  500 519  537  285 ] ;
    node_num = zeros( size(node) ) ;
    sigz = zeros( size( node ) ) ;
    for n=1:length(node)
        for ie=1:element_number
            index = find( gElement(ie,1:20) == node(n) ) ;
            if ~isempty( index )
                sigz( n ) = sigz( n ) + gElementStress( ie, index, 3 ) ;
                node_num(n) = node_num(n) + 1; 
            end
        end
    end
    sigz = sigz./node_num ;
    y = gNode( node, 2 ) ;
    sigz1 = 16e6*(0-gNode(node(1),3))/pi * y ;
    figure ;
    plot( sigz1, y, '-', sigz, y, 'o' ) ;
    xlabel( 'sigma-z' ) ;
    ylabel( 'y' ) ;
    
     % 读取剪应力沿截面宽度的变化
    node = [285 366 365 364 363 362 361 360 137] ;
    node_num = zeros( size(node) ) ;
    taoyz = zeros( size( node ) ) ;
    for n=1:length(node)
        for ie=1:element_number
            index = find( gElement(ie,1:20) == node(n) ) ;
            if ~isempty( index )
                taoyz( n ) = taoyz( n ) + gElementStress( ie, index, 4 ) ;
                node_num(n) = node_num(n) + 1; 
            end
        end
    end
    taoyz = taoyz./node_num ;
    x = gNode( node, 1 ) ;
    mu = gMaterial( 1, 2 ) ;
    taoyz1 = -(3+2*mu)/8/(1+mu)*16e6/pi*(1-(1-2*mu)/(3+2*mu)*x.^2) ;
    figure ;
    plot( x, taoyz1, '-', x, taoyz, 'o' ) ;
    xlabel( 'x' ) ;
    ylabel( 'tao-xz' ) ;
    set( gca, 'ylim', [min(taoyz), 0]*1.2 ) ;
    
    fprintf( '\n\n\n          有限元与解析解的应力比较(MPa)\n' ) ;
    fprintf( '==========================================================\n' ) ;
    fprintf( '       正应力(sigma-x)                剪应力(tao-yz)\n' ) ;
    fprintf( ' 位置(y)  有限元   解析解       位置(x)   有限元    解析解\n' ) ;
    fprintf( '----------------------------------------------------------\n' ) ;
    for n=1:length(y)
        fprintf( '  %6.3f   %5.2f   %5.2f         %6.3f    %5.2f    %5.2f\n', ...
                 y(n), sigz(n)/1e6, sigz1(n)/1e6, ...
                 x(n), taoyz(n)/1e6, taoyz1(n)/1e6 ) ;
    end
    fprintf( '==========================================================\n' ) ;
return 