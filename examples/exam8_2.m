function exam8_2
% 本程序为第八章的第二个算例，采用平面梁单元计算两铰抛物线拱的在初始条件下
%  自由振动，并对时程曲线结果进行FFT变换，求得的频率可与exam8_1.m的结果进行
%  比较，以验证本程序的可靠性
%      输入参数： 无
%      输出结果： 位移的时程曲线及其频谱特性图 

    PlaneFrameModel ;             % 定义有限元模型
    SolveModel ;                  % 求解有限元模型
    SaveResults('exam8_2.mat') ;  % 保存计算结果
    DisplayResults ;              % 显示计算结果
return ;

function PlaneFrameModel
%  定义平面杆系的有限元模型
%  输入参数：
%      无
%  返回值：
%      无
%  说明：
%      该函数定义平面杆系的有限元模型数据：
%        gNode -------- 节点定义
%        gElement ----- 单元定义
%        gMaterial ---- 材料定义，包括弹性模量，梁的截面积和梁的抗弯惯性矩
%        gBC1 --------- 约束条件
%        gDeltaT ------ 时间步长
%        gTimeEnd ----- 计算结束时刻
%        gDisp -------- 位移时程响应
%        gVelo -------- 速度时程响应
%        gAcce -------- 加速度时程响应

    global gNode gElement gMaterial gBC1 gDeltaT gTimeEnd gDisp gVelo gAcce

    % 给定抛物线拱的几何特征
    L = 60 ;               %  计算跨径(m)     
    f = 7.5 ;              %  计算矢高(m)
    
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
     
    gDeltaT = 0.01 ;
    gTimeEnd = 4096*gDeltaT  ;    % 计算时间为载荷通过所需时间的两倍
    timestep = floor(gTimeEnd/gDeltaT) ;

    % 定义位移，速度和加速度
    gDisp = zeros( (n+1)*3, timestep ) ;
    gVelo = zeros( (n+1)*3, timestep ) ;
    gAcce = zeros( (n+1)*3, timestep ) ;
    
    % 初始条件
    gDisp(:,1) = zeros( (n+1)*3, 1 ) ;
    gVelo(:,1) = ones( (n+1)*3, 1 ) ;
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
%        2. 用Newmark法计算时程响应

    global gNode gElement gMaterial gBC1 gK gM gDeltaT gTimeEnd gDisp gVelo gAcce

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

    % step3. 计算时程响应(Newmark法)
    % step3.1 初始计算
    gama = 0.5 ;
    beta = 0.25 ;
    C = zeros( size( gK ) ) ;
    [N,N] = size( gK ) ;
    alpha0 = 1/beta/gDeltaT^2 ;
    alpha1 = gama/beta/gDeltaT ;
    alpha2 = 1/beta/gDeltaT ;
    alpha3 = 1/2/beta - 1 ;
    alpha4 = gama/beta - 1 ;
    alpha5 = gDeltaT/2*(gama/beta-2) ;
    alpha6 = gDeltaT*(1-gama) ;
    alpha7 = gama*gDeltaT ;
    K1 = gK + alpha0*gM + alpha1*C;
    timestep = floor(gTimeEnd/gDeltaT) ;
    
    % step3.2 对K1进行边界条件处理
    [bc1_number,dummy] = size( gBC1 ) ;
    K1im = zeros(N,bc1_number) ;
    for ibc=1:1:bc1_number
        n = gBC1(ibc, 1 ) ;
        d = gBC1(ibc, 2 ) ;
        m = (n-1)*3 + d ;
        K1im(:,ibc) = K1(:,m) ; 
        K1(:,m) = zeros( node_number*3, 1 ) ;
        K1(m,:) = zeros( 1, node_number*3 ) ;
        K1(m,m) = 1.0 ;
    end
    [KL,KU] = lu(K1) ;   % 进行三角分解，节省后面的求解时间
    
    % step3.3 计算初始加速度
    gAcce(:,1) = gM\(-gK*gDisp(:,1)-C*gVelo(:,1)) ;
    
    % step3.4 对每一个时间步计算
    for i=2:1:timestep
        if mod(i,100) == 0
            fprintf( '当前时间步：%d\n', i ) ;
        end
        f1 = gM*(alpha0*gDisp(:,i-1)+alpha2*gVelo(:,i-1)+alpha3*gAcce(:,i-1)) ...
                  + C*(alpha1*gDisp(:,i-1)+alpha4*gVelo(:,i-1)+alpha5*gAcce(:,i-1)) ;
        % 对f1进行边界条件处理
        [bc1_number,dummy] = size( gBC1 ) ;
        for ibc=1:1:bc1_number
            n = gBC1(ibc, 1 ) ;
            d = gBC1(ibc, 2 ) ;
            m = (n-1)*3 + d ;
            f1 = f1 - gBC1(ibc,3) * K1im(:,ibc) ;
            f1(m) = gBC1(ibc,3) ;
        end
        y = KL\f1 ;
        gDisp(:,i) = KU\y ;
        gAcce(:,i) = alpha0*(gDisp(:,i)-gDisp(:,i-1)) - alpha2*gVelo(:,i-1) - alpha3*gAcce(:,i-1) ;
        gVelo(:,i) = gVelo(:,i-1) + alpha6*gAcce(:,i-1) + alpha7*gAcce(:,i) ;
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

function SaveResults( file_out )
%  保存计算结果
%  输入参数：
%     无
%  返回值：
%     无
    global gNode gElement gMaterial gBC1 gDeltaT gTimeEnd gLoad gLoadVelo gDisp gVelo gAcce
    save( file_out, 'gNode', 'gElement', 'gMaterial', 'gBC1',  ...
          'gDeltaT', 'gTimeEnd', 'gLoad', 'gLoadVelo', 'gDisp', 'gVelo', 'gAcce' ) ;
return

function DisplayResults
%  显示计算结果
%  输入参数：
%     无
%  返回值：
%     无

    global gNode gElement gMaterial gBC1 gDisp gVelo gAcce gDeltaT gTimeEnd

    % 绘制时程曲线
    [node_number,dummy] = size(gNode) ;
    t = 0:gDeltaT:gTimeEnd-gDeltaT ;
    d = gDisp((floor(node_number/4)*3)+2,:) ;
    subplot(2,1,1) ;
    plot( t, d ) ;
    title( 'L/4处挠度时程曲线' ) ;
    xlabel( '时间(s)') ;
    ylabel( '挠度(m)' ) ;
    
    % 对时程曲线进行FFT变换，获取频谱特性
    fd = fft( d ) ;
    df = 1/gTimeEnd ;
    f = (0:length(d)-1)*df ;
    subplot(2,1,2);
    plot(f,abs(fd)) ;
    set(gca,'xlim',[2,10]) ;
    title( 'L/4处挠度的频谱图' ) ;
    xlabel( '频率(Hz)') ;
    ylabel( '幅值' ) ;
    
    % 标注频率峰值
    fifi1 = diff(abs(fd));
    n = length(fifi1) ;
    d1 = fifi1(1:n-1);
    d2 = fifi1(2:n) ;
    indmax = find( d1.*d2<0 & d1>0 )+1;
    for i=1:length(indmax)
        if f(indmax(i)) > 10 
            break ;
        end
        text( f(indmax(i)+2), abs(fd(indmax(i)))*0.9, sprintf('f=%.3f',f(indmax(i))));
    end
return
