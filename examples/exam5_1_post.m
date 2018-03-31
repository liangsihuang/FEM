function [x,sy,syi]=exam5_1_post( file_mat, y, yi )
%  本程序是第5章算例1的后处理程序
%     主要是绘制单向拉伸杆中沿横截面的应力分布图
%  输入参数
%     y -------- 是一个 1×n 的数组，包含 n 个 y 坐标
%     yi ------- 是一个 1×k 的数组，包含 k 个 y 坐标
%  返回值
%     x -------- 在横截面上m个离散点
%     sy ------- 在n个横截面上对应于x的应力值
%     syi ------ 对应于yi的应力值

    global gNode gElement gMaterial gBC1 gDF gNF gK gDelta gNodeStress gElementStress
    
    if nargin < 1
        file_mat = 'exam5_1.mat' ;
        y = [3,3.5,4,4.5] ;
        yi = 0:0.25:4.5 ;
    elseif nargin < 2
        y = [3,3.5,4,4.5] ;
        yi = 0:0.25:4.5 ;
    elseif nargin < 3
        yi = 0:0.25:4.5 ;
    end
    
    load( file_mat ) ;
    
    if min(y) < 0 | max(y) > 5
        disp( 'y值超出范围，请在[0,5]区间内取值' ) ;
    end

    figure ;
    axis off ;
    axis equal ;
    str = [] ;
    sy = [] ;
    line( [-1,1,1,-1,-1], [0, 0, 5, 5, 0] ) ;
    for j=1:length(y)
        x = 0:0.05:1 ;
        z = 0.0 ;
        S = zeros(1, length(x) ) ;
        for i=1:length(x)
            disp( sprintf( 'y=%f, x=%f\n', x, y ) ) ;
            str = FindStress( [x(i),y(j),z] ) ;
            S(i) = str( 2 ) ;
       end
       line([-1,1],[y(j),y(j)]) ;
       f = 0.2/100e6 ;
       line(x,y(j)+f*S); 
       line(-x,y(j)+f*S);
       sy = [sy; S] ;
    end
    title( '单向拉伸杆正应力沿横截面的分布' ) ;
    
    syi = zeros(size(yi));
    for k=1:length(yi)
        disp( sprintf( 'k=%d\n', k ) ) ;
        str = FindStress( [0, yi(k), 0] ) ;
        syi(k) = str(2) ;
    end
    figure ;
    plot( yi, syi, '-o' ) ;
    title( '单向拉伸杆正应力的分布' ) ;
    xlabel( 'y' ) ;
    ylabel( '拉应力' ) ;

    figure ;
    hold ;
    color = 'bgrc' ;
    for i=1:length(y)
        plot( x, sy(i,:),color(i) ) ;
    end
    legend( 'y=3.0','y=3.5','y=4.0','y=4.5' ) ;
    title( '离集中力不同距离处，正应力沿横截面的分布' ) ;
    xlabel( 'x' ) ;
    ylabel( '正应力' ) ;
    hold off ;
return

function str = FindStress( pos )
%  获取某点的6个应力分量
%  输入参数
%     pos -------- 是一个 1×3 的数组，是所求点的3个坐标
%  返回值
%     str -------- 该点的6个应力分量( sx, sy, sz, tyz, txz, txy )
    global gNode gElement gMaterial gBC1 gDF gNF gK gDelta gNodeStress gElementStress
    
    % 确定点所在的单元
    [element_number,dummy] = size( gElement ) ;
    for ie=1:element_number
        if IsInElement( pos, ie )
            break ;
        end
    end
    
    % 利用体积坐标插值计算应力
    x = gNode( gElement( ie, 1:4 ), 1 )' ;
    y = gNode( gElement( ie, 1:4 ), 2 )' ;
    z = gNode( gElement( ie, 1:4 ), 3 )' ;
    V = abs( TetrahedronVol( x, y, z ) ) ;
    Vp = abs( TetrahedronVol( [x(1:3),pos(1)], [y(1:3),pos(2)], [z(1:3),pos(3)] ) ) ;
    Vi = abs( TetrahedronVol( [x(2),x(4),x(3),pos(1)], [y(2),y(4),y(3),pos(2)],[z(2),z(4),z(3),pos(3)] ) );
    Vj = abs( TetrahedronVol( [x(3),x(4),x(1),pos(1)], [y(3),y(4),y(1),pos(2)],[z(3),z(4),z(1),pos(3)] ) );
    Vm = abs( TetrahedronVol( [x(1),x(4),x(2),pos(1)], [y(1),y(4),y(2),pos(2)],[z(1),z(4),z(2),pos(3)] ) );
    Si = Vi/V ;
    Sj = Vj/V ;
    Sm = Vm/V ;
    Sp = Vp/V ;
    str =  Si.*gNodeStress(gElement(ie,1),:) + Sj.*gNodeStress(gElement(ie,2),:) ...
         + Sm.*gNodeStress(gElement(ie,3),:) + Sp.*gNodeStress(gElement(ie,4),:) ;
return

function b = IsInElement( pos, ie )
%  确定空间点是否在一个单元内
%  输入参数
%     pos -------- 是一个 1×3 的数组，是空间点的3个坐标
%     ie --------- 单元号
%  返回值
%     b ---------- 可以是下列值之一
%                   0  ==>  单元外
%                   1  ==>  单元内
    global gNode gElement

    % 获取单元的4个结点坐标，
    i = gElement( ie, 1 ) ;
    j = gElement( ie, 2 ) ;
    m = gElement( ie, 3 ) ;
    p = gElement( ie, 4 ) ;
    
    x = gNode( gElement(ie,1:4), 1 ) ;
    y = gNode( gElement(ie,1:4), 2 ) ;
    z = gNode( gElement(ie,1:4), 3 ) ;
    
    xmin = min(x) ;
    xmax = max(x) ;
    if pos(1) > xmax | pos(1) < xmin
        b = 0 ;
        return ;
    end
    
    ymin = min(y) ;
    ymax = max(y) ;
    if pos(2) > ymax | pos(2) < ymin
        b = 0 ;
        return ;
    end
    
    zmin = min(z) ;
    zmax = max(z) ;
    if pos(3) > zmax | pos(3) < zmin
        b = 0 ;
        return ;
    end
    
    % 如果四个结点的顺序不符合右手系，则交换它们的顺序
    v = TetrahedronVol( x, y, z ) ;
    if v < 0
        k = j ;
        j = m ;
        m = k ;
    end
    
    % 如果该空间点在单元内，则4个小四面体的体积都不小于零
    delta = [ i, j, m ;
              j, p, m ;
              m, p, i ;
              i, p, j ] ;
    for k=1:4
        x = gNode( delta(k,:), 1 ) ; 
        y = gNode( delta(k,:), 2 ) ;
        z = gNode( delta(k,:), 3 ) ;
        x = [x',pos(1)] ;
        y = [y',pos(2)] ;
        z = [z',pos(3)] ;
        v = TetrahedronVol( x, y, z ) ;
        if v < 0
            b = 0 ;
            return
        end
    end
    b = 1 ;
return

function v = TetrahedronVol( x, y, z )
%  计算一个四面体的体积
%  输入参数
%     x -------- 四面体4个顶点的x坐标
%     y -------- 四面体4个顶点的y坐标
%     z -------- 四面体4个顶点的z坐标
%  返回值
%     v -------- 四面体的体积，如果是负值，表示这四个顶点的顺序没有按右手系定义

    v = det( [ 1 x(1) y(1) z(1)
               1 x(2) y(2) z(2)
               1 x(3) y(3) z(3)
               1 x(4) y(4) z(4) ] ) ;
    v = v / 6 ;
return