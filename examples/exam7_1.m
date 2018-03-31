function exam7_1
% 本程序为第七章的第一个算例，采用薄板矩形单元计算矩形板的弯曲
%  输入参数： 
%      无

    % 定义板的尺寸
    a = 1 ;     % 板x方向宽度
    b = 1 ;     % 板y方向长度
    h = 0.01 ;  % 板z方向厚度
    
    % 定义网格密度     
    m = 16 ;     % x方向单元数目 
    n = 16 ;     % y方向单元数目

    % 定义材料性质
    E = 2.1e11 ;   % 弹性模量
    mu = 0.3 ;     % 泊松比
    
    % 荷载参数
    load = '均布力' ;   % 也可以是 集中力
    %load = '集中力' ;
    p = 1e6 ;             % 载荷大小
    
    % 边界条件
    %bc = '简支' ;  % 也可以是 固支    
    bc = '固支' ;
    
    % 输出文件名称
    file_out = 'exam7_1.mat' ;
    
    % 检查输出文件是否存在
    file_prompt = 0 ;
    if exist( file_out ) ~= 0 & file_prompt == 1
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

    FemModel(a,b,h,m,n,E,mu,load,p,bc) ;       % 定义有限元模型
    SolveModel ;                     % 求解有限元模型
    SaveResults( file_out ) ;        % 保存计算结果
    DisplayResults(a,b,h,m,n,E,mu,load,p,bc) ; % 把计算结果显示出来 
return

function FemModel(a,b,h,m,n,E,mu,load,p,bc)
%  定义有限元模型
%  输入参数：
%    a  ---------   板x方向宽度
%    b  ---------   板y方向长度
%    h  ---------   板z方向厚度
%    m  ---------   x方向单元数目 
%    n  ---------   y方向单元数目
%    E  ---------   弹性模量
%    mu ---------   泊松比
%    load -------   载荷类型
%    p  ---------   载荷大小
%    bc ---------   边界条件
%  返回值：
%      无

    global gNode gElement gMaterial gBC1 gDF gNF gK

    % 计算结点坐标
    dx = a / m / 2 ;   % 取板的1/4进行分析
    dy = b / n / 2 ;
    gNode = zeros( (m+1)*(n+1), 2 ) ;
    for i=1:n+1
        for j=1:m+1
            gNode( (i-1)*(m+1)+j, : ) = [dx*(j-1),dy*(i-1)] ;
        end
    end
    
    % 定义单元
    gElement = zeros( n*n, 5 ) ;
    for i=1:n
        for j=1:m
            gElement( (i-1)*m+j, 1:4) = [ (i-1)*(m+1)+j, ...
                                          (i-1)*(m+1)+j+1, ...
                                          i*(m+1)+j+1,...
                                          i*(m+1)+j ] ;
        end
    end
    gElement( :, 5 ) = 1 ;
    
    % 定义材料
    gMaterial = [ E, mu, h] ;  
    
    % 确定边界条件
    if bc == '简支'
        gBC1 = zeros( (n+1)*2+m*2+1, 3 ) ;
        for i=1:m+1
            gBC1( i, : ) = [i, 2, 0] ;                    % x=0的边界上绕x轴的转角等于零
        end
        for i=1:n+1
            gBC1( (m+1)+i, : ) = [(i-1)*(m+1)+1, 3, 0] ;  % y=0的边界上绕y轴的转角等于零
        end
        for i=1:n+1
            gBC1( m+n+2+i, : ) = [i*(m+1), 1, 0 ] ;       % x=a/2的边界上挠度为零
        end
        for i=1:m
            gBC1( (n+1)*2+m+1+i, : ) = [n*(m+1)+i, 1, 0] ;    % y=b/2的边界上挠度为零
        end
    elseif bc == '固支'
        gBC1 = zeros( m*4+(n+1)*3+n, 3 ) ;
        for i=1:m
            gBC1( i, : ) = [i, 2, 0] ;                    % x=0的边界上绕x轴的转角等于零
        end
        for i=1:n
            gBC1( m+i, : ) = [(i-1)*(m+1)+1, 3, 0] ;      % y=0的边界上绕y轴的转角等于零
        end
        for i=1:n+1
            gBC1( m+n+i, : ) = [i*(m+1), 1, 0 ] ;          % x=a/2的边界上挠度为零
        end
        for i=1:n+1
            gBC1( m+n+n+1+i, : ) = [i*(m+1), 2, 0 ] ;      % x=a/2的边界上绕x轴的转角为零
        end
        for i=1:n+1
            gBC1( m+n+(n+1)*2+i, : ) = [i*(m+1), 3, 0 ] ;  % x=a/2的边界上绕y轴的转角为零
        end
        for i=1:m
            gBC1( m+n+(n+1)*3+i, : ) = [n*(m+1)+i, 1, 0] ; % y=b/2的边界上挠度为零
        end
        for i=1:m
            gBC1( m*2+n+(n+1)*3+i, : ) = [n*(m+1)+i, 2, 0] ; % y=b/2的边界上绕x轴的转角为零
        end
        for i=1:m
            gBC1( m*3+n+(n+1)*3+i, : ) = [n*(m+1)+i, 3, 0] ; % y=b/2的边界上绕y轴的转角为零
        end
    end

    % 确定载荷
    if load == '均布力'
        gNF = [] ;
        gDF = gElement ;
        gDF(:,5) = p ;
    elseif load == '集中力'
        gDF = [] ;
        gNF = [1, 1, p] ;
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

    global gNode gElement gMaterial gBC1 gDF gNF gK gDelta gElementMoment gNodeMoment

    % step1. 定义整体刚度矩阵和节点力向量
    [node_number,dummy] = size( gNode ) ;
    gK = sparse( node_number * 3, node_number * 3 ) ;
    f = sparse( node_number * 3, 1) ;    

    % step2. 计算单元刚度矩阵，并集成到整体刚度矩阵中
    [element_number,dummy] = size( gElement ) ;
    for ie=1:1:element_number
        disp( sprintf(  '计算单元刚度矩阵并集成，当前单元: %d', ie  ) ) ;
        k = StiffnessMatrix( ie ) ;
        AssembleStiffnessMatrix( ie, k ) ;
    end
    
    % step3. 计算分布压力的等效节点力，并集成到整体节点力向量中
    [df_number, dummy] = size( gDF ) ;
    for idf = 1:1:df_number
        edf = EquivalentDistPressure( gDF(idf,1:4), gDF(idf,5) ) ;
        node = gDF( idf, 1:4 ) ;
        f( node*3-2 ) = f( node*3-2 ) + edf( 1:3:10 ) ;
        f( node*3-1 ) = f( node*3-1 ) + edf( 2:3:11 ) ;
        f( node*3 )   = f( node*3 )   + edf( 3:3:12 ) ;
    end
  
    % step4. 计算分布压力的等效节点力，并集成到整体节点力向量中
    [nf_number,dummy] = size( gNF ) ;
    for inf = 1:nf_number
        node = gNF( inf, 1 ) ;
        dim = gNF( inf, 2 ) ;
        f( (node-1)*3+dim ) = f( (node-1)*3+dim ) + gNF( inf, 3 ) ;
    end
    % step5. 处理第一类约束条件，修改刚度矩阵和节点力向量。采用乘大数法
    [bc1_number,dummy] = size( gBC1 ) ;
    for ibc=1:1:bc1_number
        n = gBC1(ibc, 1 ) ;
        d = gBC1(ibc, 2 ) ;
        m = (n-1)*3 + d ;
        f(m) = gBC1(ibc, 3)* gK(m,m) * 1e10 ;
        gK(m,m) = gK(m,m) * 1e10 ;
    end
    
    % step5. 求解方程组，得到节点位移向量
    gDelta = gK \ f ;
    
    % step6. 计算单元弯矩
    gElementMoment = zeros( element_number, 4, 3 ) ;
    for ie=1:element_number 
        disp( sprintf(  '计算单元弯矩，当前单元: %d', ie  ) ) ;
        em = ElementMoment( ie ) ;
        gElementMoment( ie, :, : ) = em ;
    end
    
    % step7. 计算节点弯矩(采用绕节点平均)
    gNodeMoment = zeros( node_number, 3 ) ;       
    nodenum = zeros( node_number,1 ) ;
    for i=1:node_number                         
        disp( sprintf(  '计算节点弯矩，当前结点: %d', i  ) ) ;
        for ie=1:element_number
            index = find(  gElement( ie, 1:4 ) == i ) ;
            if ~isempty(index)
                nodenum(i) = nodenum(i) + 1 ;
                gNodeMoment(i,1) = gNodeMoment(i,1) + gElementMoment( ie, index, 1 ) ;
                gNodeMoment(i,2) = gNodeMoment(i,2) + gElementMoment( ie, index, 2 ) ;
                gNodeMoment(i,3) = gNodeMoment(i,3) + gElementMoment( ie, index, 3 ) ;
            end
        end
    end
    gNodeMoment(:,1) = gNodeMoment(:,1) ./ nodenum ;
    gNodeMoment(:,2) = gNodeMoment(:,2) ./ nodenum ;
    gNodeMoment(:,3) = gNodeMoment(:,3) ./ nodenum ;
return

function k = StiffnessMatrix( ie )   
%  计算平面应变等参数单元的刚度矩阵
%  输入参数：
%     ie -- 单元号
%  返回值：
%     k  -- 单元刚度矩阵

    global gNode gElement gMaterial
    k = zeros( 12, 12 ) ;
    E = gMaterial( gElement( ie, 5 ), 1 ) ;
    mu = gMaterial( gElement( ie, 5 ), 2 ) ;
    h = gMaterial( gElement( ie, 5 ), 3 ) ;
    D = E * h^3 / 12 / ( 1 - mu^2 ) ;
    x = gNode( gElement( ie, : ), 1 ) ;
    y = gNode( gElement( ie, : ), 2 ) ;
    a = abs( x(3) - x(1) ) / 2 ;
    b = abs( y(3) - y(1) ) / 2 ;
    ab2 = a^2/b^2 ;
    ba2 = b^2/a^2 ;
    H = D/60/a/b ;
    xi = [ -1, 1, 1, -1 ] ;
    eta = [ -1, -1, 1, 1 ] ;
    for i=1:4
        for j=1:4
            xi0 = xi(i)*xi(j) ;
            eta0 = eta(i)*eta(j) ;
            a11 = 3*H*(15*(ba2*xi0+ab2*eta0) ...
                   +(14-4*mu+5*ba2+5*ab2)*xi0*eta0) ;
            a12 = -3*H*b*((2+3*mu+5*ab2)*xi0*eta(i) ...
                   + 15*ab2*eta(i)+5*mu*xi0*eta(j)) ;
            a13 = 3*H*a*((2+3*mu+5*ba2)*xi(i)*eta0 ...
                   + 15*ba2*xi(i)+5*mu*xi(j)*eta0) ;
            a21 = -3*H*b*((2+3*mu+5*ab2)*xi0*eta(j) ...
                   + 15*ab2*eta(j)+5*mu*xi0*eta(i)) ;
            a22 = H*b^2*(2*(1-mu)*xi0*(3+5*eta0)+...
                   + 5*ab2*(3+xi0)*(3+eta0)) ;
            a23 = -15*H*mu*a*b*(xi(i)+xi(j))*(eta(i)+eta(j)) ;
            a31 = 3*H*a*((2+3*mu+5*ba2)*xi(j)*eta0 ...
                   + 15*ba2*xi(j)+5*mu*xi(i)*eta0) ;
            a32 = -15*H*mu*a*b*(xi(i)+xi(j))*(eta(i)+eta(j)) ;
            a33 = H*a^2*(2*(1-mu)*eta0*(3+5*xi0)+...
                   + 5*ba2*(3+xi0)*(3+eta0)) ;
            k( i*3-2:i*3, j*3-2:j*3 ) = [ a11 a12 a13;
                                          a21 a22 a23;
                                          a31 a32 a33 ] ;
        end
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

function edf = EquivalentDistPressure( node, press )
%   计算分布压力的等效节点力
%   输入参数:
%      node ----------  结点号
%      press ---------  跟结点号对应的压力值
%   返回值:
%      edf -----------  等效节点力向量
    global gNode
    
    edf = zeros( 12, 1 ) ;
    x = gNode( node, 1 ) ;
    y = gNode( node, 2 ) ;
    a = abs( x(3)-x(1) ) / 2 ;
    b = abs( y(3)-y(1) ) / 2 ;
    xi = [-1; 1; 1; -1] ;
    eta = [-1;-1; 1; 1] ;
    for i=1:4
        edf(i*3-2:i*3) = [press*a*b; -press*a*b^2/3*eta(i); press*a^2*b/3*xi(i)] ;
    end
return

function em = ElementMoment( ie )
%  计算单元的结点弯矩
%  输入参数:
%      ie ----- 单元号
%  返回值:
%      em ----- 单元结点弯矩

    global gElement gDelta gMaterial

    E = gMaterial( gElement( ie, 5 ),  1 ) ;
    mu = gMaterial( gElement( ie, 5 ), 2 ) ;
    h = gMaterial( gElement( ie, 5 ),  3 ) ;
    em = zeros( 4, 3 ) ;
    node = gElement( ie, 1:4 ) ;
    delta = zeros( 12, 1 ) ;
    delta(1:3:10) = gDelta( (node-1)*3+1 ) ;
    delta(2:3:11) = gDelta( (node-1)*3+2 ) ;
    delta(3:3:12) = gDelta( (node-1)*3+3 ) ;
    xi  = [ -1, 1, 1, -1 ] ;
    eta = [ -1, -1, 1, 1 ] ;
    D = E/(1-mu^2)*[1  mu  0
                    mu  1  0
                    0   0  (1-mu)/2 ] ;
    for i=1:4
        B = MatrixB( ie, xi(i), eta(i) ) ;
        em( i, : ) = transpose(h^3/12*D*B*delta) ;
    end
return

function B = MatrixB( ie, x, e )
%  计算单元的应变矩阵B
%  输入参数:
%      ie ----- 单元号
%  返回值:
%      B ------ 单元应变矩阵
    
    global gNode gElement
    B = zeros( 3, 12 ) ;
    xy = gNode( gElement(ie, 1:4), : ) ;
    a = abs(xy(3,1)-xy(1,1))/2 ;
    b = abs(xy(3,2)-xy(1,2))/2 ;
    xi  = [ -1, 1, 1, -1 ] ;
    eta = [ -1, -1, 1, 1 ] ;
    for i=1:4
        x0 = xi(i)*x ;
        e0 = eta(i)*e ;
        B(:,(i-1)*3+1:(i-1)*3+3) = [ 3*b/a*x0*(1+e0)                       0                      b*xi(i)*(1+3*x0)*(1+e0)
                                     3*a/b*e0*(1+x0)                -a*eta(i)*(1+x0)*(1+3*e0)           0
                                     xi(i)*eta(i)*(3*x^2+3*e^2-4)   -b*xi(i)*(3*e^2+2*e0-1)       a*eta(i)*(3*x^2+2*x0-1) ] ;
    end
    B = B / a / b / 4 ;
return

function SaveResults( file_out )
%  保存计算结果
%  输入参数：
%     无
%  返回值：
%     无
    global gNode gElement gMaterial gBC1 gDF gDelta gNF gK
    save( file_out, 'gNode', 'gElement', 'gMaterial', 'gBC1',  ...
          'gDF', 'gDelta', 'gNF', 'gK' ) ;
return

function DisplayResults(a,b,h,m,n,E,mu,load,p,bc)
%  显示计算结果，并与解析解结果比较
%  输入参数：
%    a  ---------   板x方向宽度
%    b  ---------   板y方向长度
%    h  ---------   板z方向厚度
%    m  ---------   x方向单元数目 
%    n  ---------   y方向单元数目
%    E  ---------   弹性模量
%    mu ---------   泊松比
%    load -------   载荷类型
%    p  ---------   载荷大小
%    bc ---------   边界条件
%  返回值：
%     无

    global gDelta gNodeMoment
    
    wmax = gDelta( 1, 1 ) ;          % 中心挠度
    D = E*h^3/12/(1-mu^2) ;          % 抗弯刚度

    mx  = gNodeMoment( 1, 1 ) ;   % 中心内力矩(Mx)
    my  = gNodeMoment( 1, 2 ) ;   % 中心内力矩(My)
    mxy = gNodeMoment( 1, 3 ) ;   % 中心内力矩(Mxy)
    
    fprintf( '板的尺寸a×b×h = %g×%g×%g\n', a, b, h ) ;
    fprintf( '弹性模量 = %g\n', E ) ;
    fprintf( '泊松比   = %g\n', mu ) ;
    fprintf( '载荷类型 = %s\n', load ) ;
    fprintf( '载荷大小 = %g\n', p ) ;
    fprintf( '边界条件 = %s\n', bc ) ;
    fprintf( '网格密度 m×n = %d×%d\n', m, n ) ;
    if load == '均布力'
        alpha = wmax * D / p / a^4 ;
        beta = mx/p/a^2 ;
        beta1 = my/p/a^2 ;
    elseif load == '集中力'
        alpha = wmax * D / p / 4 / a^2 ;
        beta = mx/p/4 ;
        beta1 = my/p/3 ;
    end
    
    fprintf( '挠度系数 alpha = %.7f\n', alpha ) ;
    fprintf( '弯矩系数 beta  = %.7f\n', beta ) ;
    fprintf( '弯矩系数 beta1 = %.7f\n', beta1 ) ;
%    fid = fopen( 'exam7_1.out', 'a' ) ;
%        fprintf( fid, '%.1f  %15.7f  %15.7f\n', b/a, alpha, beta ) ;
%    fclose( fid ) ;
return