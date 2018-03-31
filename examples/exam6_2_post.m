    if exist( 'exam6_2.mat' ) == 0
        disp( sprintf( '�����ļ� %s ������', file_in ) )
        disp( sprintf( '������exam6_2.m��������exam6_2.mat' ) ) ;
        disp( sprintf( '������ֹ' ) )
        return ;
    end
    load exam6_2.mat
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
    subplot( 1, 2, 1 ) ;
    plot( sigz1/1e6, y, '-', sigz/1e6, y, 'o' ) ;
    xlabel( 'sigma-z (MPa) ' ) ;
    ylabel( 'y' ) ;
    legend( '������', '����Ԫ' ) ;
    title( '������Ӧ�������ߵķֲ�ͼ' ) ;
    
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
    subplot( 1, 2, 2 ) ;
    plot( x, taoyz1/1e6, '-', x, taoyz/1e6, 'o' ) ;
    xlabel( 'x' ) ;
    ylabel( 'tao-xz (MPa)' ) ;
    set( gca, 'ylim', [min(taoyz/1e6), 0]*1.2 ) ;
    legend( '������', '����Ԫ' ) ;
    title( '����������Ӧ��������ķֲ�ͼ' ) ;
    set(gcf, 'Position', [100,100,600,400]);
    
    fprintf( '\n\n\n          ����Ԫ��������Ӧ���Ƚ�(MPa)\n' ) ;
    fprintf( '==========================================================\n' ) ;
    fprintf( '       ��Ӧ��(sigma-x)                ��Ӧ��(tao-yz)\n' ) ;
    fprintf( ' λ��(y)  ����Ԫ   ������       λ��(x)   ����Ԫ    ������\n' ) ;
    fprintf( '----------------------------------------------------------\n' ) ;
    for n=1:length(y)
        fprintf( '  %6.3f   %5.2f   %5.2f         %6.3f    %5.2f    %5.2f\n', ...
                 y(n), sigz(n)/1e6, sigz1(n)/1e6, ...
                 x(n), taoyz(n)/1e6, taoyz1(n)/1e6 ) ;
    end
    fprintf( '==========================================================\n' ) ;
