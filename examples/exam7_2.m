function exam7_2
% 本程序为第七章的第二个算例，采用薄板三角形单元计算斜板的弯曲
%  输入参数： 
%       无

    %   斜板尺寸示意图
    %                        a 
    %            |-----------------------|                  a :  长度(x)  
    %                                                       b :  宽度(y)  
    %       -|-   \-----------------------\                 h :  厚度
    %        |     \                       \   |           fi :  斜角(单位: 度)
    %     b  |      \                       \fi|
    %        |       \                       \-|
    %       -|-       \-----------------------\|    
    
    % 定义板的尺寸 
    a = 1.92 ;     
    b = 1 ;        
    h = 0.01 ;     
    fi = 30 ;      

    % 定义网格密度     
    m = 16 ;     % 长度划分的数目 
    n = 16 ;     % 宽度划分的数目

    % 定义材料性质
    E = 2.1e11 ;   % 弹性模量
    mu = 0.2 ;     % 泊松比
    
    % 定义均布荷载大小
    p = 1.0e6 ;

    % 输出文件名称
    file_out = 'exam7_2.mat' ;
    
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

    FemModel(a,b,fi,h,m,n,E,mu,p) ;         % 定义有限元模型
    SolveModel ;                            % 求解有限元模型
    SaveResults( file_out ) ;               % 保存计算结果
    DisplayResults(a,b,fi,h,m,n,E,mu,p) ;   % 把计算结果显示出来 
return

function FemModel(a,b,fi,h,m,n,E,mu,p)
%  定义有限元模型
%  输入参数：
%    a  ---------   板长度
%    b  ---------   板宽度
%    fi ---------   斜角
%    h  ---------   板厚度
%    m  ---------   长度方向单元数目 
%    n  ---------   宽度方向单元数目
%    E  ---------   弹性模量
%    mu ---------   泊松比
%    p ----------   均布载荷大小
%  返回值：
%      无

    global gNode gElement gMaterial gBC1 gDF
    
    dx =  a / m ;
    dy =  b / n ;
    tanfi = tan(fi*pi/180) ;
    
    gNode = zeros( (n+1)*(m+1), 2 ) ;
    for i=1:n+1
        for j=1:m+1
            gNode( (i-1)*(m+1)+j, : ) = [(b-(i-1)*dy)*tanfi+dx*(j-1),dy*(i-1)] ;
        end
    end
    
    gElement = zeros( n*m*2, 4 ) ;
    for i=1:n
        for j=1:m
            gElement( (i-1)*m*2+(j-1)*2+1, 1:3) = [ (i-1)*(m+1)+j, ...
                                                    i*(m+1)+j+1, ...
                                                    i*(m+1)+j ] ;
            gElement( (i-1)*m*2+(j-1)*2+2, 1:3) = [ (i-1)*(m+1)+j, ...
                                                    (i-1)*(m+1)+j+1, ...
                                                    i*(m+1)+j+1 ] ;
        end
    end
    gElement( :, 4 ) = 1 ;
    
    gMaterial = [ E, mu, h] ;  % 材料特性
    
    gBC1 = zeros( (n+1)*2, 3 ) ;
    for i=1:n+1
        gBC1( i, : ) = [ (i-1)*(m+1)+1, 1, 0 ] ;
    end
    for i=1:n+1
        gBC1( (n+1)+i, : ) = [i*(m+1), 1, 0] ;
    end

    gDF = gElement ;
    gDF( :, 4 ) = p ;
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

    global gNode gElement gMaterial gBC1 gDF gK gDelta gElementMoment gNodeMoment

    % step1. 定义整体刚度矩阵和节点力向量
    [node_number,dummy] = size( gNode ) ;
    gK = sparse( node_number * 3, node_number * 3 ) ;
    f = sparse( node_number * 3, 1) ;    

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
        edf = EquivalentDistPressure( gDF(idf,1:3), gDF(idf,4) ) ;
        node = gDF( idf, 1:3) ;
        f( node*3-2 ) = f( node*3-2 ) + edf( 1:3:7 ) ;
        f( node*3-1 ) = f( node*3-1 ) + edf( 2:3:8 ) ;
        f( node*3 )   = f( node*3 )   + edf( 3:3:9 ) ;
    end
  
    % step4. 处理第一类约束条件，修改刚度矩阵和节点力向量。采用乘大数法
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
    
    % step6. 计算单元的每个结点的弯矩
    gElementMoment = zeros( element_number, 3, 3 ) ;
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
            index = find(  gElement( ie, 1:3 ) == i ) ;
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
%  计算薄板三角形单元的刚度矩阵
%  输入参数：
%     ie -- 单元号
%  返回值：
%     k  -- 单元刚度矩阵
    global gNode gElement
    
    x = gNode( gElement( ie, 1:3 ), 1 ) ;
    y = gNode( gElement( ie, 1:3 ), 2 ) ;
    i=1; j=2; m=3;
    ai = x(j)*y(m) - x(m)*y(j) ;
    aj = x(m)*y(i) - x(i)*y(m) ;
    am = x(i)*y(j) - x(j)*y(i) ;
    delta2 = ai+aj+am ;
    
    k = zeros( 9, 9 ) ;
    L = [1/2 1/2   0;
           0 1/2 1/2;
         1/2   0 1/2] ;
    W = [1/6 1/6 1/6] ;
    D = MatrixD( ie ) ;
    for i=1:length(W)
         B = MatrixB( ie, L(i,:) ) ;
         k = k + W(i)*transpose(B)*D*B*delta2 ;
    end
return

function B = MatrixB( ie, L )
%  计算薄板三角形单元的应变矩阵
%  输入参数：
%     ie -- 单元号
%  返回值：
%     B  -- 单元应变矩阵
    global gNode gElement
    
    x = gNode( gElement( ie, 1:3 ), 1 ) ;
    y = gNode( gElement( ie, 1:3 ), 2 ) ;
    i=1; j=2; m=3;
    ai = x(j)*y(m) - x(m)*y(j) ;
    aj = x(m)*y(i) - x(i)*y(m) ;
    am = x(i)*y(j) - x(j)*y(i) ;
    bi = y(j) - y(m) ;
    bj = y(m) - y(i) ;
    bm = y(i) - y(j) ;
    ci = -(x(j)-x(m)) ;
    cj = -(x(m)-x(i)) ;
    cm = -(x(i)-x(j)) ;
    delta2 = ai+aj+am ;
    T = [    bi^2          bj^2             2*bi*bj
             ci^2          cj^2             2*ci*cj
          2*bi*ci       2*bj*cj     2*(bi*cj+bj*ci) ] ;
    Li = L(1); Lj = L(2); Lm = L(3) ;
    Nii = [ -4*Lj-12*Li+6, ...
            -6*bj*Li+2*bj-3*bj*Lj-bm*Lj, ...
            2*cj-6*cj*Li-3*cj*Lj-cm*Lj, ...
            -4*Lj, ...
            (-bm+bi)*Lj, ...
            -(cm-ci)*Lj, ...
            -6+12*Li+8*Lj, ...
            bi*Lj-6*bj*Li+4*bj-3*bj*Lj, ...
            ci*Lj+4*cj-6*cj*Li-3*cj*Lj ] ;
    Njj = [ -4*Li, ...
            -(bj-bm)*Li, ...
            (-cj+cm)*Li, ...
            6-4*Li-12*Lj, ...
            bm*Li+6*bi*Lj-2*bi+3*bi*Li, ...
            cm*Li+6*ci*Lj-2*ci+3*ci*Li, ...
            8*Li-6+12*Lj, ...
            6*bi*Lj-4*bi+3*bi*Li-bj*Li, ...
            6*ci*Lj-4*ci+3*ci*Li-cj*Li] ;
    Nij = [ -4*Li-4*Lj+2, ...
            -3*bj*Li-bm*Li+1/2*bj-bj*Lj-1/2*bm+bm*Lj, ...
            -3*cj*Li-cm*Li+1/2*cj-cj*Lj-1/2*cm+cm*Lj, ...
            -4*Li-4*Lj+2, ...
            bm*Lj+3*bi*Lj+1/2*bm-bm*Li-1/2*bi+bi*Li, ...
            cm*Lj+3*ci*Lj+1/2*cm-cm*Li-1/2*ci+ci*Li, ...
            -4+8*Li+8*Lj, ...
            3*bi*Lj-3/2*bi+bi*Li-3*bj*Li+3/2*bj-bj*Lj, ...
            3*ci*Lj-3/2*ci+ci*Li-3*cj*Li+3/2*cj-cj*Lj] ;
    B = -1/delta2^2*T*[Nii; Njj; Nij] ;
return

function D = MatrixD( ie )
%  计算薄板三角形单元的弹性矩阵
%  输入参数：
%     ie -- 单元号
%  返回值：
%     D  -- 单元弹性矩阵
    global gElement gMaterial
    E = gMaterial( gElement( ie, 4 ), 1 ) ;
    mu = gMaterial( gElement( ie, 4 ), 2 ) ;
    h = gMaterial( gElement( ie, 4 ), 3 ) ;
    D = E*h^3/12/(1-mu^2)*[  1        mu           0
                            mu         1           0
                             0         0    (1-mu)/2 ] ;
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
    edf = zeros( 9, 1 ) ;
    x = gNode( node, 1 ) ;
    y = gNode( node, 2 ) ;
    i=1; j=2; m=3;
    ai = x(j)*y(m) - x(m)*y(j) ;
    aj = x(m)*y(i) - x(i)*y(m) ;
    am = x(i)*y(j) - x(j)*y(i) ;
    bi = y(j) - y(m) ;
    bj = y(m) - y(i) ;
    bm = y(i) - y(j) ;
    ci = -(x(j)-x(m)) ;
    cj = -(x(m)-x(i)) ;
    cm = -(x(i)-x(j)) ;
    delta = abs(ai+aj+am)/2 ;
    edf = delta*press/24*[ 8; bj-bm; cj-cm;
                           8; bm-bi; cm-ci;
                           8; bi-bj; ci-cj] ;
return

function em = ElementMoment( ie )
%  计算单元的结点弯矩
%  输入参数:
%      ie ----- 单元号
%  返回值:
%      em ----- 单元结点弯矩

    global gElement gDelta gMaterial

    E = gMaterial( gElement( ie, 4 ),  1 ) ;
    mu = gMaterial( gElement( ie, 4 ), 2 ) ;
    h = gMaterial( gElement( ie, 4 ),  3 ) ;
    em = zeros( 3, 3 ) ;
    node = gElement( ie, 1:3 ) ;
    delta = zeros( 9, 1 ) ;
    delta(1:3:7) = gDelta( (node-1)*3+1 ) ;
    delta(2:3:8) = gDelta( (node-1)*3+2 ) ;
    delta(3:3:9) = gDelta( (node-1)*3+3 ) ;
    L  = [ 1, 0, 0
           0, 1, 0
           0, 0, 1] ;
    D = E/(1-mu^2)*[1  mu  0
                    mu  1  0
                    0   0  (1-mu)/2 ] ;
    for i=1:3
        B = MatrixB( ie, L(i,:) ) ;
        em( i, : ) = transpose(h^3/12*D*B*delta) ;
    end
return
    
function SaveResults( file_out )
%  保存计算结果
%  输入参数：
%     无
%  返回值：
%     无
    global gNode gElement gMaterial gBC1 gDF gDelta gElementMoment gNodeMoment
    save( file_out, 'gNode', 'gElement', 'gMaterial', 'gBC1',  ...
          'gDF', 'gDelta', 'gElementMoment', 'gNodeMoment' ) ;
return

function DisplayResults(a,b,fi,h,m,n,E,mu,p)
%  显示计算结果，并与解析解结果比较
%  输入参数：
%    a  ---------   板长度
%    b  ---------   板宽度
%    fi ---------   斜角
%    h  ---------   板厚度
%    m  ---------   长度方向单元数目 
%    n  ---------   宽度方向单元数目
%    E  ---------   弹性模量
%    mu ---------   泊松比
%    p ----------   均布载荷大小
%  返回值：
%     无
    global gMaterial gNode gDelta gElementMoment gNodeMoment
    [node_number, dummy] = size( gNode ) ;
    for i=1:node_number
        fprintf( '%6d  %10.5f  %10.5f  %18.8e  %18.8e  %18.8e\n', i, gNode(i,1), gNode(i,2),...
               gDelta( i*3-2 ), gDelta(i*3-1), gDelta(i*3) ) ;
    end
    
    % 计算板中心的挠度和内力矩
    node0 = (n/2+1)*(m+1)+(m/2+1) ;
    w0 = gDelta( node0*3-2 ) ;
    Mx0 = gNodeMoment( node0, 1 ) ;
    My0 = gNodeMoment( node0, 2 ) ;
    Mxy0 = gNodeMoment( node0, 3 ) ;
    M0 = Mx0*cos(fi*pi/180)^2+My0*sin(fi*pi/180)^2+Mxy0*sin(2*fi*pi/180) ;
    
    % 计算自由边上最大挠度和最大弯矩(方法是利用三次样条加密，然后取最大值)
    x = gNode(1:m+1,1) ;
    w1 = full(gDelta(1:3:(m+1)*3-2)) ;
    M1 = gNodeMoment(1:m+1,1) ;
    xx = x(1):0.01:x(length(x)) ;
    ww = spline(x,w1,xx) ;
    MM = spline(x,M1,xx) ;
    w1max = max(ww) ;
    M1max = max(MM) ;
    
    % 计算挠度和内力矩系数，并显示
    D = E*h^3/12/(1-mu^2) ;
    alpha0 = w0*D/p/b^4 ;
    alpha1 = w1max*D/p/b^4 ;
    beta0 = M0/p/b^2 ;
    beta1 = M1max/p/b^2 ;
    fprintf( 'w0 = %g    alpha0 = %.6f\n', w0, alpha0 ) ;
    fprintf( 'w1 = %g    alpha1 = %.6f\n', w1max, alpha1 ) ;
    fprintf( 'M0 = %g    beta0 = %.6f\n', M0, beta0 ) ;
    fprintf( 'M1 = %g    beta1 = %.6f\n', M1max, beta1 ) ;
return