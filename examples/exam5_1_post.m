function [x,sy,syi]=exam5_1_post( file_mat, y, yi )
%  �������ǵ�5������1�ĺ������
%     ��Ҫ�ǻ��Ƶ�����������غ�����Ӧ���ֲ�ͼ
%  �������
%     y -------- ��һ�� 1��n �����飬���� n �� y ����
%     yi ------- ��һ�� 1��k �����飬���� k �� y ����
%  ����ֵ
%     x -------- �ں������m����ɢ��
%     sy ------- ��n��������϶�Ӧ��x��Ӧ��ֵ
%     syi ------ ��Ӧ��yi��Ӧ��ֵ

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
        disp( 'yֵ������Χ������[0,5]������ȡֵ' ) ;
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
    title( '�����������Ӧ���غ����ķֲ�' ) ;
    
    syi = zeros(size(yi));
    for k=1:length(yi)
        disp( sprintf( 'k=%d\n', k ) ) ;
        str = FindStress( [0, yi(k), 0] ) ;
        syi(k) = str(2) ;
    end
    figure ;
    plot( yi, syi, '-o' ) ;
    title( '�����������Ӧ���ķֲ�' ) ;
    xlabel( 'y' ) ;
    ylabel( '��Ӧ��' ) ;

    figure ;
    hold ;
    color = 'bgrc' ;
    for i=1:length(y)
        plot( x, sy(i,:),color(i) ) ;
    end
    legend( 'y=3.0','y=3.5','y=4.0','y=4.5' ) ;
    title( '�뼯������ͬ���봦����Ӧ���غ����ķֲ�' ) ;
    xlabel( 'x' ) ;
    ylabel( '��Ӧ��' ) ;
    hold off ;
return

function str = FindStress( pos )
%  ��ȡĳ���6��Ӧ������
%  �������
%     pos -------- ��һ�� 1��3 �����飬��������3������
%  ����ֵ
%     str -------- �õ��6��Ӧ������( sx, sy, sz, tyz, txz, txy )
    global gNode gElement gMaterial gBC1 gDF gNF gK gDelta gNodeStress gElementStress
    
    % ȷ�������ڵĵ�Ԫ
    [element_number,dummy] = size( gElement ) ;
    for ie=1:element_number
        if IsInElement( pos, ie )
            break ;
        end
    end
    
    % ������������ֵ����Ӧ��
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
%  ȷ���ռ���Ƿ���һ����Ԫ��
%  �������
%     pos -------- ��һ�� 1��3 �����飬�ǿռ���3������
%     ie --------- ��Ԫ��
%  ����ֵ
%     b ---------- ����������ֵ֮һ
%                   0  ==>  ��Ԫ��
%                   1  ==>  ��Ԫ��
    global gNode gElement

    % ��ȡ��Ԫ��4��������꣬
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
    
    % ����ĸ�����˳�򲻷�������ϵ���򽻻����ǵ�˳��
    v = TetrahedronVol( x, y, z ) ;
    if v < 0
        k = j ;
        j = m ;
        m = k ;
    end
    
    % ����ÿռ���ڵ�Ԫ�ڣ���4��С��������������С����
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
%  ����һ������������
%  �������
%     x -------- ������4�������x����
%     y -------- ������4�������y����
%     z -------- ������4�������z����
%  ����ֵ
%     v -------- ����������������Ǹ�ֵ����ʾ���ĸ������˳��û�а�����ϵ����

    v = det( [ 1 x(1) y(1) z(1)
               1 x(2) y(2) z(2)
               1 x(3) y(3) z(3)
               1 x(4) y(4) z(4) ] ) ;
    v = v / 6 ;
return