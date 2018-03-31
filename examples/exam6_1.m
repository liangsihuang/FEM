function exam6_1
% 本程序为第六章的第一个算例，采用4～8结点四边形等参元计算受内压的旋转厚壁圆筒
%  输入参数： 
%       无

    % 输出文件名称
    file_out = 'exam6_1.mat' ;
    
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

    FemModel ;                  % 定义有限元模型
    SolveModel ;                % 求解有限元模型
    SaveResults( file_out ) ;   % 保存计算结果
    DisplayResults ;            % 把计算结果显示出来 
return

function FemModel
%  定义有限元模型
%  输入参数：
%      无
%  返回值：
%      无

%  全局变量及名称
%      gNode ------ 节点坐标
%      gElement --- 单元定义
%      gMaterial -- 材料性质
%      gBC1 ------- 第一类约束条件
%      gBC3 ------- 斜约束
%      gDF -------- 分布力
    global gNode gElement gMaterial gBC1 gBC3 gDF
    
    r1 = 0.05 ;              %  内半径
    r2 = 0.1 ;               %  外半径
    theta = 10.0*pi/180.0 ;  %  圆心角（弧度）
    
    % 节点坐标
    node_number = 28 ;
    gNode = zeros( node_number, 2 ) ;
    j = 0 ;
    for i = 1:5:26     % 结点1,6,11,16,21,26的坐标
        r = r1 + (r2-r1)/5 * j ;
        gNode(i,1) = r*cos(theta) ;
        gNode(i,2) = r*sin(theta) ;
        j = j + 1 ;
    end
    j = 0 ;
    for i = 4:5:24     % 结点4,9,14,19,24的坐标
        r = r1+(r2-r1)/10+(r2-r1)/5*j ;
        gNode(i,1) = r*cos(theta) ;
        gNode(i,2) = r*sin(theta) ;
        j = j+1 ;
    end
    j = 0 ;
    for i=2:5:27       % 结点2,7,12,17,22,27的坐标
        r = r1+(r2-r1)/5*j ;
        gNode(i,1) = r*cos(theta/2) ;
        gNode(i,2) = r*sin(theta/2) ;
        j = j+1 ;
    end
    j = 0 ;
    for i=3:5:28      % 结点3,8,13,18,23,28的坐标
        r = r1+(r2-r1)/5*j ;
        gNode(i,1) = r ;
        gNode(i,2) = 0 ;
        j = j+1 ;
    end
    j = 0 ;
    for i=5:5:25      % 结点5,10,15,20,25的坐标
        r = r1+(r2-r1)/10+(r2-r1)/5*j ;
        gNode(i,1) = r ;
        gNode(i,2) = 0 ;
        j = j+1 ;
    end
    
    % 单元定义
    element_number = 5 ;
    gElement = zeros( element_number, 9 ) ;
    gElement(1,:)=[3  8  6  1   5  7  4  2  1] ;
    for i=2:5
        gElement(i,1:8) = gElement(i-1,1:8)+5 ;
        gElement(i,9) = 1 ;
    end
    
    % 材料性质
    gMaterial = [2.0e11 0.3 7800] ;
    
    % 第一类约束条件
    bc1_number = 11;
    gBC1 = zeros( bc1_number, 3 ) ;
    j = 0 ;
    for i=3:5:28      %  结点3,8,13,18,23,28的边界条件(y方向约束)
        j = j + 1 ;
        gBC1(j,1) = i ;  % 结点号
        gBC1(j,2) = 2 ;  % 自由度号
        gBC1(j,3) = 0 ;  % 约束值
    end
    j = 6 ;
    for i=5:5:25     %  结点5,10,15,20,25的边界条件(y方向约束)
        j = j + 1;
        gBC1(j,1) = i ;  % 结点号
        gBC1(j,2) = 2 ;  % 自由度号
        gBC1(j,3) = 0 ;  % 约束值
    end
    
    % 斜约束条件
    bc3_number = 17 ;
    gBC3 = zeros( bc3_number, 4 ) ;
    j = 0 ;
    for i=1:5:26      %  结点1,6,11,16,21,26的边界条件(y'方向约束)
        j = j + 1 ;
        gBC3(j,1) = i ;     % 结点号
        gBC3(j,2) = 2 ;     % 自由度号
        gBC3(j,3) = 0 ;     % 约束值
        gBC3(j,4) = theta ; % 倾斜角度
    end
    j = 6 ;
    for i=4:5:24     %  结点4,9,14,19,24的边界条件(y'方向约束)
        j = j + 1;
        gBC3(j,1) = i ;     % 结点号
        gBC3(j,2) = 2 ;     % 自由度号
        gBC3(j,3) = 0 ;     % 约束值
        gBC3(j,4) = theta ; % 倾斜角度
    end
    j = 11 ;
    for i=2:5:27 
        j = j+1 ;
        gBC3(j,1) = i ;       % 结点号
        gBC3(j,2) = 2 ;       % 自由度号
        gBC3(j,3) = 0 ;       % 约束值
        gBC3(j,4) = theta/2 ; % 倾斜角度
    end
    
    % 分布压力
    df_number = 1 ;
    gDF = [3 1 2 120e6 120e6 120e6] ;
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

    global gNode gElement gMaterial gBC1 gBC3 gDF gK gDelta1 gDelta gElementStress

    % step1. 定义整体刚度矩阵和节点力向量
    [node_number,dummy] = size( gNode ) ;
    gK = sparse( node_number * 2, node_number * 2 ) ;
    f = sparse( node_number * 2, 2 ) ;    

    % step2. 计算单元刚度矩阵，并集成到整体刚度矩阵中
    [element_number,dummy] = size( gElement ) ;
    for ie=1:1:element_number
        disp( sprintf(  '计算单元刚度矩阵并集成，当前单元: %d', ie  ) ) ;
        k = StiffnessMatrix( ie ) ;
        AssembleStiffnessMatrix( ie, k ) ;
    end
    
    % step3. 计算单元离心力的等效结点力，并集成到整体结点力向量中
    for ie=1:1:element_number
        evf = EquivalentVolumeForce( ie ) ;
        node = gElement( ie, 1:8 ) ;
        f( (node-1)*2+1,2 ) = f( (node-1)*2+1,2 ) + evf(1:2:15) ;
        f( (node-1)*2+2,2 ) = f( (node-1)*2+2,2 ) + evf(2:2:16) ;
    end

    % step4. 计算分布压力的等效节点力，并集成到整体节点力向量中
    [df_number, dummy] = size( gDF ) ;
    for idf = 1:1:df_number
        enf = EquivalentDistPressure( gDF(idf,1:3), gDF(idf,4:6) ) ;
        i = gDF(idf, 1) ;
        j = gDF(idf, 2) ;
        m = gDF(idf, 3) ;
        f( (i-1)*2+1:(i-1)*2+2, 1 ) = f( (i-1)*2+1:(i-1)*2+2, 1 ) + enf( 1:2 ) ;
        f( (j-1)*2+1:(j-1)*2+2, 1 ) = f( (j-1)*2+1:(j-1)*2+2, 1 ) + enf( 3:4 ) ;
        f( (m-1)*2+1:(m-1)*2+2, 1 ) = f( (m-1)*2+1:(m-1)*2+2, 1 ) + enf( 5:6 ) ;
    end
  
    % step5. 处理斜约束边界条件
    [bc3_number,dummy] = size(gBC3) ;
    for ibc=1:1:bc3_number
        n = gBC3( ibc, 1 ) ;
        theta = gBC3( ibc, 4 ) ;
        c = cos(theta) ;
        s = sin(theta) ;
        T = [ c s
             -s c ] ;
        gK((n-1)*2+1:(n-1)*2+2,:) = T*gK((n-1)*2+1:(n-1)*2+2,:) ;
        gK(:,(n-1)*2+1:(n-1)*2+2) = gK(:,(n-1)*2+1:(n-1)*2+2)*transpose(T) ;
        f((n-1)*2+1:(n-1)*2+2,:) = T*f((n-1)*2+1:(n-1)*2+2,:) ;
        gBC1 = [gBC1; gBC3(ibc,1:3)] ;
    end
    
    % step6. 处理第一类约束条件，修改刚度矩阵和节点力向量。采用乘大数法
    [bc1_number,dummy] = size( gBC1 ) ;
    for ibc=1:1:bc1_number
        n = gBC1(ibc, 1 ) ;
        d = gBC1(ibc, 2 ) ;
        m = (n-1)*2 + d ;
        f(m,:) = gBC1(ibc, 3)* gK(m,m) * 1e15 ;
        gK(m,m) = gK(m,m) * 1e15 ;
    end
    
    % step7. 求解方程组，得到节点位移向量
    gDelta1 = gK \ f ;
    
    % step8. 把有斜约束的结点位移转换到整体坐标
    gDelta = gDelta1 ;
    for ibc=1:bc3_number
        n = gBC3( ibc, 1 ) ;
        theta = gBC3( ibc, 4 ) ;
        c = cos(theta) ;
        s = sin(theta) ;
        T = [ c s
             -s c ] ;
        gDelta( (n-1)*2+1:n*2,:) = transpose(T)*gDelta( (n-1)*2+1:n*2,: ) ;
    end
    
    % step9. 计算每个单元的结点应力
    gElementStress = zeros( element_number, 8, 3, 2 ) ;
    delta = zeros( 16, 2 ) ;
    for ie = 1:element_number
        xi  = [ -1   1   1  -1   0   1   0  -1 ] ;
        eta = [ -1  -1   1   1  -1   0   1   0 ] ;
        for n=1:8
            B = MatrixB( ie, xi(n), eta(n) ) ;
            D = MatrixD( ie, 1 ) ;
            delta(1:2:15, :) = gDelta( gElement(ie,1:8)*2-1, : ) ;
            delta(2:2:16, :) = gDelta( gElement(ie,1:8)*2, : ) ;
            sigma = D*B*delta ;
            gElementStress( ie, n, :, : ) = sigma ;
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

    k = zeros( 16, 16 ) ;
    D = MatrixD( ie, 1 ) ;
    % 3 x 3  高斯积分点和权系数
    %x = [-0.774596669241483,              0.0,0.774596669241483] ;   
    %w = [ 0.555555555555556,0.888888888888889,0.555555555555556] ;  
    % 2 x 2  高斯积分点和权系数
    x = [-0.577350269189626, 0.577350269189626] ;                    
    w = [1, 1] ;                                                     
    for i=1:1:length(x)
        for j=1:1:length(x)
            B = MatrixB( ie, x(i), x(j) ) ;
            J = Jacobi( ie, x(i), x(j) ) ;
            k = k + w(i)*w(j)*transpose(B)*D*B*det(J) ;   
        end
    end
return

function D = MatrixD( ie, opt )
%  计算单元的弹性矩阵D
%  输入参数：
%     ie --------- 单元号
%    opt --------- 问题选项
%                   1 ==>  平面应力
%                   2 ==>  平面应变
%  返回值：
%     D  --------- 弹性矩阵D

    global gElement gMaterial 
    E  = gMaterial( gElement(ie, 9), 1 ) ;   % 弹性模量
    mu = gMaterial( gElement(ie, 9), 2 ) ;   % 泊松比 
    if opt == 1   % 平面应力的弹性常数
        A1 = mu ;                                
        A2 = (1-mu)/2 ;                          
        A3  = E/(1-mu^2) ;                      
    else          % 平面应变的弹性常数
       A1 = mu/(1-mu) ;                        
        A2 = (1-2*mu)/2/(1-mu) ;                
        A3 = E*(1-mu)/(1+mu)/(1-2*mu) ;         
    end
    D = A3* [  1  A1   0
              A1   1   0
               0   0  A2] ;
return

function B = MatrixB( ie, xi, eta )
%  计算单元的应变矩阵B
%  输入参数：
%     ie --------- 单元号
%     xi,eta ----- 局部坐标  
%  返回值：
%     B  --------- 在局部坐标处的应变矩阵B

    [N_x,N_y] = N_xy( ie, xi, eta ); 
    B = zeros( 3, 16 ) ;
    for i=1:1:8
        B(1:3,(2*i-1):2*i) = [ N_x(i)        0      
                                    0   N_y(i)
                               N_y(i),  N_x(i)]; 
    end
return

function [N_x, N_y] = N_xy( ie, xi, eta )
%  计算形函数对整体坐标的导数
%  输入参数：
%     ie --------- 单元号
%     xi,eta ----- 局部坐标  
%  返回值：
%     N_x  ------- 在局部坐标处的形函数对x坐标的导数
%     N_y  ------- 在局部坐标处的形函数对y坐标的导数

    J = Jacobi( ie, xi, eta ) ;
    [N_xi,N_eta] = N_xieta( ie, xi, eta ) ;
    A=inv(J)*[N_xi;N_eta] ;
    N_x = A(1,:) ;
    N_y = A(2,:) ;
return

function [N_xi, N_eta] = N_xieta( ie, xi, eta )
%  计算形函数对局部坐标的导数
%  输入参数：
%     ie --------- 单元号
%     xi,eta ----- 局部坐标  
%  返回值：
%     N_xi   ------- 在局部坐标处的形函数对xi坐标的导数
%     N_eta  ------- 在局部坐标处的形函数对eta坐标的导数

    x = [ -1, 1, 1, -1 ] ;
    e = [ -1, -1, 1, 1 ] ;
    N_xi  = zeros( 1, 8 ) ;
    N_eta = zeros( 1, 8 ) ;

    N_xi( 5 )  = xi*(eta-1) ;
    N_eta( 5 ) = 0.5*(xi^2-1) ;
    N_xi( 6 )  = 0.5*(1-eta^2) ;
    N_eta( 6 ) = -eta*(xi+1) ;
    N_xi( 7 )  = -xi*(eta+1) ;
    N_eta( 7 ) = 0.5*(1-xi^2) ;
    N_xi( 8 )  = 0.5*(eta^2-1) ;
    N_eta( 8 ) = eta*(xi-1) ;

    N_xi(1)  = x(1)*(1+e(1)*eta)/4 - 0.5*( N_xi(5)  + N_xi(8) );
    N_eta(1) = e(1)*(1+x(1)*xi)/4  - 0.5*( N_eta(5) + N_eta(8) ) ;
    N_xi(2)  = x(2)*(1+e(2)*eta)/4 - 0.5*( N_xi(5)  + N_xi(6) );
    N_eta(2) = e(2)*(1+x(2)*xi)/4  - 0.5*( N_eta(5) + N_eta(6) ) ;
    N_xi(3)  = x(3)*(1+e(3)*eta)/4 - 0.5*( N_xi(6)  + N_xi(7) );
    N_eta(3) = e(3)*(1+x(3)*xi)/4  - 0.5*( N_eta(6) + N_eta(7) ) ;
    N_xi(4)  = x(4)*(1+e(4)*eta)/4 - 0.5*( N_xi(7)  + N_xi(8) );
    N_eta(4) = e(4)*(1+x(4)*xi)/4  - 0.5*( N_eta(7) + N_eta(8) ) ;
return

function J = Jacobi( ie, xi, eta )
%  计算雅克比矩阵
%  输入参数：
%     ie --------- 单元号
%     xi,eta ----- 局部坐标  
%  返回值：
%     J   ------- 在局部坐标(xi,eta)处的雅克比矩阵
    global gNode gElement
    x = gNode(gElement(ie,1:8),1) ;
    y = gNode(gElement(ie,1:8),2) ;
    [N_xi,N_eta] = N_xieta( ie, xi, eta ) ;
    x_xi  = N_xi  * x ;
    x_eta = N_eta * x ;
    y_xi  = N_xi  * y ;
    y_eta = N_eta * y ;
    J = [ x_xi, y_xi; x_eta, y_eta ];
return

function AssembleStiffnessMatrix( ie, k )
%  把单元刚度矩阵集成到整体刚度矩阵
%  输入参数:
%      ie  --- 单元号
%      k   --- 单元刚度矩阵
%  返回值:
%      无
    global gElement gK
    for i=1:1:8
        for j=1:1:8
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

function edf = EquivalentDistPressure( node, pressure )
%   计算分布压力的等效节点力
%   输入参数:
%      node ----------  结点号
%      pressure ------  跟结点号对应的压力值
%   返回值:
%      edf -----------  等效节点力向量
    global gNode
    x = gNode( node, 1 ) ;
    y = gNode( node, 2 ) ;
    xi = x(1) ; xj = x(2) ; xk = x(3) ;
    yi = y(1) ; yj = y(2) ; yk = y(3) ;
    sigi = pressure(1) ;
    sigj = pressure(2) ;
    sigk = pressure(3) ;
    
    X1 = -10*xi-2*xj+12*xk ;
    Y1 = -10*yi-2*yj+12*yk ;
    X2 = xi - xj ;
    Y2 = yi - yj ;
    X3 = -6*xi-2*xj+8*xk ;
    Y3 = -6*yi-2*yj+8*yk ;
    X4 = 2*xi+10*xj-12*xk ;
    Y4 = 2*yi+10*yj-12*yk ;
    X5 = 2*xi+6*xj-8*xk ;
    Y5 = 2*yi+6*yj-8*yk ;
    X6 = -16*xi+16*xj ;
    Y6 = -16*yi+16*yj ;
    
    XY = [  Y1  Y2  Y3
           -X1 -X2 -X3
            Y2  Y4  Y5
           -X2 -X4 -X5
            Y3  Y5  Y6
           -X3 -X5 -X6 ] ;
    sig = [sigi; sigj; sigk ] ;
    edf = 1/30*XY*sig ;
return

function evf = EquivalentVolumeForce( ie )
%   计算离心力的等效节点力
%   输入参数:
%      ie ----------  单元号
%   返回值:
%      evf ---------  等效节点力向量
    global gNode gElement gMaterial
    
    evf = zeros( 16, 1 ) ;
    omega = 2094.0 ;   % 旋转角速度
    ro = gMaterial( gElement( ie, 9 ), 3 ) ;
    x = gNode( gElement(ie,1:8), 1 ) ;
    y = gNode( gElement(ie,1:8), 2 ) ;
    xi = [-0.774596669241483,0.0,              0.774596669241483] ;
    w  = [ 0.555555555555556,0.888888888888889,0.555555555555556] ;
    for i=1:length(xi)
        for j=1:length(xi)
            J = Jacobi( ie, xi(i), xi(j) ) ;
            detJ = det(J); 
            N = ShapeFunction( xi(i), xi(j) ) ;
            evf(1:2:15) = evf(1:2:15) + N'*N*x*detJ*w(i)*w(j) ;
            evf(2:2:16) = evf(2:2:16) + N'*N*y*detJ*w(i)*w(j) ;
        end
    end
    evf = evf * ro * omega^2 ;
return

function N = ShapeFunction( xi, eta )
%   计算形函数的值
%   输入参数:
%      ie ----------  单元号
%      xi, eta -----  单元内局部坐标
%   返回值:
%      N -----------  形函数的值

    N5 = ( eta - 1 ) * ( xi^2 - 1 ) / 2 ;
    N6 = ( xi + 1 ) * ( 1 - eta^2 ) / 2 ;
    N7 = ( eta + 1 ) * ( 1 - xi^2 ) / 2 ;
    N8 = ( xi - 1 ) * ( eta^2 - 1 ) / 2 ;
    N1 = ( 1 - xi ) * ( 1 - eta ) / 4 - 0.5 * ( N8 + N5 ) ;
    N2 = ( 1 + xi ) * ( 1 - eta ) / 4 - 0.5 * ( N5 + N6 ) ;
    N3 = ( 1 + xi ) * ( 1 + eta ) / 4 - 0.5 * ( N6 + N7 ) ;
    N4 = ( 1 - xi ) * ( 1 + eta ) / 4 - 0.5 * ( N7 + N8 ) ;
    N = [ N1 N2 N3 N4 N5 N6 N7 N8 ]; 
return

function SaveResults( file_out ) 
%  保存计算结果
%  输入参数：
%     无
%  返回值：
%     无
    global gNode gElement gMaterial gBC1 gBC3 gDF gK gDelta gDelta1 gElementStress
    save( file_out, 'gNode', 'gElement', 'gMaterial', 'gBC1', 'gBC3', ...
          'gDF', 'gDelta', 'gDelta1', 'gElementStress' ) ;
return

function DisplayResults
%  显示计算结果，并与解析解结果比较
%  输入参数：
%     无
%  返回值：
%     无

%  说明:
%    1。受内压的轴对称圆筒平面应力问题的解析解是(参考文献[1])
%           径向位移:  ur = -C1/E*(1+mu)/r + 2*C2/E*(1-mu)*r
%                      C1 = -a^2*b^2*p/(b^2-a^2)
%                      C2 = p*a^2/2/(b^2-a^2)
%           径向应力:  sr = a^2*p/(b^2-a^2)*(1-b^2/r^2)
%           周向应力:  st = a^2*p/(b^2-a^2)*(1+b^2/r^2)
%    2。以匀角速度w旋转的圆筒平面应力问题的解析解(参考文献[2])
%           径向位移:  ur = 1/E*((1-mu)*C3*r - (1+mu)*C4/r - (1-mu^2)/8*ro*w^2*r^3 ) 
%                      C3 = (3+mu)/8*ro*w^2*(b^2+a^2)
%                      C4 = -(3+mu)/8*ro*w^2*a^2*b^2
%           径向应力:  sr = (3+mu)/8*ro*w^2*(b^2+a^2-a^2*b^2/r^2-r^2)
%           周向应力:  st = (3+mu)/8*ro*w^2*(b^2+a^2+a^2*b^2/r^2-(1+3*mu)/(3+mu)*r^2)
%     
%   参数说明
%        E  ----- 弹性模量
%       mu  ----- 泊松比
%        a  ----- 内半径
%        b  ----- 外半径
%        r  ----- 径向坐标
%        p  ----- 内压
%       ro  ----- 密度
%        w  ----- 旋转角速度(rad/s)
%  
%   参考文献
%     [1] 钱伟长 叶开元, 1956, 弹性力学, pp.230-231
%     [2] Timoshenko S. P., Goodier J. N., 1970 Theory of elasticity. New York: McGraw-Hill, pp.80-83

    global gNode gElement gMaterial gDelta gDelta1 gElementStress

    E  = gMaterial( 1, 1 ) ;
    mu = gMaterial( 1, 2 ) ;
    ro = gMaterial( 1, 3 ) ;
    a = 0.05 ;
    b = 0.10 ;
    p = 120e6 ;
    w = 2094 ;
    r = 0.05:0.005:0.1 ;

    %  计算受内压圆筒的位移和应力的解析解
    C1 = -a^2*b^2*p/(b^2-a^2) ;
    C2 = p*a^2/2/(b^2-a^2) ;
    ur1 = -C1/E*(1+mu)./r + 2*C2/E*(1-mu)*r ;
    sr1 = a^2*p/(b^2-a^2)*(1-b^2./r.^2) ;
    st1 = a^2*p/(b^2-a^2)*(1+b^2./r.^2) ;
    
    %  计算旋转圆盘的位移和应力的解析解
    C3 = (3+mu)/8*ro*w^2*(b^2+a^2) ;
    C4 = -(3+mu)/8*ro*w^2*a^2*b^2 ;
    ur2 = 1/E*((1-mu)*C3*r - (1+mu)*C4./r - (1-mu^2)/8*ro*w^2*r.^3 ) ;
    sr2 = (3+mu)/8*ro*w^2*(b^2+a^2-a^2*b^2./r.^2-r.^2) ;
    st2 = (3+mu)/8*ro*w^2*(b^2+a^2+a^2*b^2./r.^2-(1+3*mu)/(3+mu)*r.^2) ;
    
    %  读取有限元的计算结果
    node = [3, 5, 8, 10, 13, 15, 18, 20, 23, 25, 28]; 
    node_num = zeros( size(node) ) ;
    ur1_fem = gDelta( node*2-1, 1 ) ;
    ur2_fem = gDelta( node*2-1, 2 ) ;
    [element_number,dummy] = size( gElement ) ;
    sr1_fem = zeros( size( sr1 ) ) ;
    st1_fem = zeros( size( st1 ) ) ;
    sr2_fem = zeros( size( sr2 ) ) ;
    st2_fem = zeros( size( st2 ) ) ;
    for n=1:length(node)
        for ie=1:element_number
            index = find( gElement( ie, 1:8 ) == node( n ) ) ;
            if ~isempty(index)    
                sr1_fem( n ) = sr1_fem( n ) + gElementStress( ie, index, 1, 1 ) ;
                st1_fem( n ) = st1_fem( n ) + gElementStress( ie, index, 2, 1 ) ;
                sr2_fem( n ) = sr2_fem( n ) + gElementStress( ie, index, 1, 2 ) ;
                st2_fem( n ) = st2_fem( n ) + gElementStress( ie, index, 2, 2 ) ;
                node_num(n) = node_num(n) + 1 ;
            end
        end
    end
    sr1_fem = sr1_fem ./ node_num ;
    st1_fem = st1_fem ./ node_num ;
    sr2_fem = sr2_fem ./ node_num ;
    st2_fem = st2_fem ./ node_num ;
    
    fprintf( '\n\n\n' ) ;
    fprintf( '                        表1 径向位移 (mm)\n' ) ;
    fprintf( '==================================================================\n' ) ;
    fprintf( '结点号  坐标         承受内压                     承受离心力\n' ) ;
    fprintf( '        (mm)    有限元解      解析解        有限元解     解析解\n' ) ;
    fprintf( '------------------------------------------------------------------\n' ) ;
    for i=1:length(node)
        fprintf( '%4d   %4.0f     %.6f     %.6f       %.6f     %.6f\n', ...
             node(i), r(i)*1000, ur1_fem(i)*1000, ur1(i)*1000, ur2_fem(i)*1000, ur2(i)*1000 ) ;
    end
    fprintf( '==================================================================\n\n\n' ) ;
    

    fprintf( '                      表2 径向应力 (MPa)\n' ) ;
    fprintf( '==================================================================\n' ) ;
    fprintf( '结点号  坐标         承受内压                     承受离心力\n' ) ;
    fprintf( '        (mm)    有限元解      解析解        有限元解     解析解\n' ) ;
    fprintf( '------------------------------------------------------------------\n' ) ;
    for i=1:length(node)
        fprintf( '%4d   %4.0f    %8.2f    %8.2f       %8.2f     %8.2f\n', ...
             node(i), r(i)*1000, sr1_fem(i)/1e6, sr1(i)/1e6, sr2_fem(i)/1e6, sr2(i)/1e6 ) ;
    end
    fprintf( '==================================================================\n\n\n' ) ;

    
    fprintf( '                      表3 周向应力 (MPa)\n' ) ;
    fprintf( '==================================================================\n' ) ;
    fprintf( '结点号  坐标         承受内压                     承受离心力\n' ) ;
    fprintf( '        (mm)    有限元解      解析解        有限元解     解析解\n' ) ;
    fprintf( '------------------------------------------------------------------\n' ) ;
    for i=1:length(node)
        fprintf( '%4d   %4.0f    %8.2f    %8.2f       %8.2f     %8.2f\n', ...
             node(i), r(i)*1000, st1_fem(i)/1e6, st1(i)/1e6, st2_fem(i)/1e6, st2(i)/1e6 ) ;
    end
    fprintf( '==================================================================\n\n\n' ) ;
return