function exam4_2_post( file_in )
%  �����ε�Ԫ�ĺ������
%  �������
%     file_in  -----  �������ļ����� exam4_2.m ����
%  ����ֵ
%     ��
    global gNode gElement gMaterial gBC1 gDelta gNodeStress gElementStress
    
    if nargin < 1 
        file_in = 'exam4_2.mat' ;
    end
    
    load( file_in ) ;
    gDelta = full( gDelta ) ;
    DisplayModel ;            %  ��ʾ��ʼģ��
    for i=1:6
        PlotStress( i ) ;     %  ��ʾӦ����ͼ    
    end
    PlotDisplacement( 1 ) ;   %  ��ʾx����λ����ͼ
    PlotDisplacement( 2 ) ;   %  ��ʾy����λ����ͼ
return

function DisplayModel
%  ��ͼ�η�ʽ��ʾ����Ԫģ��
%  ���������
%     ��
%  ����ֵ��
%     ��
    global gNode gElement gMaterial gBC1 
    
    figure ;
    axis equal ;
    axis off ;
    set( gcf, 'NumberTitle', 'off' ) ;
    set( gcf, 'Name', '����Ԫģ��' ) ;
    
    % ���ݲ�ͬ�Ĳ��ϣ���ʾ��Ԫ��ɫ
    [element_number, dummy] = size( gElement ) ;
    material_color = [ 'r','g','b','c','m','y','w','k'] ;
    for i=1:element_number
        x1 = gNode( gElement( i, 1 ), 1 ) ;
        x2 = gNode( gElement( i, 2 ), 1 ) ;
        x3 = gNode( gElement( i, 3 ), 1 ) ;
        y1 = gNode( gElement( i, 1 ), 2 ) ;
        y2 = gNode( gElement( i, 2 ), 2 ) ;
        y3 = gNode( gElement( i, 3 ), 2 ) ;
        color_index = mod( gElement( i, 4 ), length( material_color ) ) ; 
        if color_index == 0 
            color_index = length( material_color ) ;
        end
        patch( [x1;x2;x3], [y1;y2;y3], material_color( color_index ) ) ;
    end
    
    % ��ʾ�߽�����
    DisplayBC( 'blue' ) ;
    
return

function DisplayBC( color )
%  ��ͼ�η�ʽ��ʾ����Ԫģ�͵ı߽�����
%  ���������
%     color  ----  �߽���������ɫ
%  ����ֵ��
%     ��
    global gNode gBC1
    
    % ȷ���߽������Ĵ�С
    xmin = min( gNode(:,1) ) ;
    xmax = max( gNode(:,1) ) ;
    factor = ( xmax - xmin ) / 25 ;
    
    [bc1_number,dummy] = size( gBC1 ) ;
    dBCSize = factor ;
    for i=1:bc1_number
        if( gBC1( i, 2 ) == 1 )  %  x����Լ��
            x0 = gNode( gBC1( i, 1 ), 1 ) ;
            y0 = gNode( gBC1( i, 1 ), 2 ) ;
            x1 = x0 - dBCSize ;
            y1 = y0 + dBCSize/2 ;
            x2 = x1 ;
            y2 = y0 - dBCSize/2 ;
            hLine = line( [x0 x1 x2 x0], [y0 y1 y2 y0] ) ;
            set( hLine, 'Color', color ) ;                
            
            xCenter = x1 - dBCSize/6 ;
            yCenter = y0 + dBCSize/4 ;
            radius = dBCSize/6 ;
            theta=0:pi/6:2*pi ;
            x = radius * cos( theta ) ;
            y = radius * sin( theta ) ;
            hLine = line( x+xCenter, y+yCenter ) ;
            set( hLine, 'Color', color ) ;
            
            hLine = line( x+xCenter, y+yCenter-dBCSize/2 ) ;
            set( hLine, 'Color', color ) ;          
            
            x0 = x0 - dBCSize - dBCSize/3 ;                   
            y0 = y0 + dBCSize/2 ;
            x1 = x0 ;
            y1 = y0 - dBCSize ;
            hLine = line( [x0, x1], [y0, y1] ) ;    
            set( hLine, 'Color', color ) ;          
            
            x = [x0 x0-dBCSize/6] ;
            y = [y0 y0-dBCSize/6] ;
            hLine = line( x, y ) ;
            set( hLine, 'Color', color ) ;
            for j=1:1:4
                hLine = line( x, y - dBCSize/4*j ); 
                set( hLine, 'Color', color ) ;
            end
        else                      %  y����Լ��
            x0 = gNode( gBC1( i, 1 ), 1 ) ;
            y0 = gNode( gBC1( i, 1 ), 2 ) ;
            x1 = x0 - dBCSize/2 ;
            y1 = y0 - dBCSize ;
            x2 = x1 + dBCSize ;
            y2 = y1 ;
            hLine = line( [x0 x1 x2 x0], [y0 y1 y2 y0] ) ; 
            set( hLine, 'Color', color ) ;      
            
            xCenter = x0 - dBCSize/4 ;
            yCenter = y1 - dBCSize/6 ;
            radius = dBCSize/6 ;
            theta=0:pi/6:2*pi ;
            x = radius * cos( theta ) ;
            y = radius * sin( theta ) ;
            hLine = line( x+xCenter, y+yCenter ) ;
            set( hLine, 'Color', color ) ; 
            
            hLine = line( x+xCenter+dBCSize/2, y+yCenter ) ;
            set( hLine, 'Color', color ) ; 
            
            hLine = line( [x1, x1+dBCSize], [y1-dBCSize/3, y1-dBCSize/3] ) ;
            set( hLine, 'Color', color ) ; 
            
            x = [x1 x1-dBCSize/6] ;
            y = [y1-dBCSize/3 y1-dBCSize/2] ;
            hLine = line( x, y ) ;
            set( hLine, 'Color', color ) ;
            for j=1:1:4
                hLine = line( x+dBCSize/4*j, y ); 
                set( hLine, 'Color', color ) ;
            end
        end
    end
return

function PlotStress( iStress )
%  ��ʾӦ����ͼ
%  ���������
%     iStress --- Ӧ������ָʾ���������������ֵ
%                 1  --  x������Ӧ��
%                 2  --  y������Ӧ��
%                 3  --  ��Ӧ��
%                 4  --  ��һ��Ӧ��
%                 5  --  �ڶ���Ӧ��
%                 6  --  ����Ӧ��
%  ����ֵ��
%     ��
    global gNode gElement gNodeStress
    
    switch iStress
    case 1
        title = ' x ������Ӧ��' ;
    case 2
        title = ' y ������Ӧ��' ;
    case 3
        title = ' ��Ӧ��' ;
    case 4
        title = ' �����Ӧ��' ;
    case 5
        title = ' ��С��Ӧ��' ;
    case 6
        title = ' ����Ӧ��' ;
    end
    
    figure ;
    axis equal ;
    axis off ;
    set( gcf, 'NumberTitle', 'off' ) ;
    set( gcf, 'Name', title ) ;
    
    stressMin = min( gNodeStress( :, iStress ) ) ;
    stressMax = max( gNodeStress( :, iStress ) ) ;
    caxis( [stressMin, stressMax] ) ;
    colormap( 'jet' ) ;

    [element_number, dummy] = size( gElement ) ;
    for ie=1:1:element_number
        x = [ gNode( gElement( ie, 1 ), 1 ) ;
              gNode( gElement( ie, 2 ), 1 ) ;
              gNode( gElement( ie, 3 ), 1 ) ] ;
        y = [ gNode( gElement( ie, 1 ), 2 ) ;
              gNode( gElement( ie, 2 ), 2 ) ;
              gNode( gElement( ie, 3 ), 2 ) ] ;
        c = [ gNodeStress( gElement( ie, 1 ), iStress ) ;
              gNodeStress( gElement( ie, 2 ), iStress ) ;
              gNodeStress( gElement( ie, 3 ), iStress ) ] ;
        set( patch( x, y, c ), 'EdgeColor', 'interp' ) ;
    end
    
    yTick = stressMin:(stressMax-stressMin)/10:stressMax ;
    Label = cell( 1, length(yTick) ); 
    for i=1:length(yTick)
        Label{i} = sprintf( '%.2fMPa', yTick(i)/1e6 ) ;
    end
    set( colorbar( 'vert' ), 'YTick', yTick, 'YTickLabelMode', 'Manual', 'YTickLabel', Label ) ;
    
    % PlotStressContour( iStress, 10, 'white' ) ;
return 

function PlotStressContour( iStress, nContour, color )
%  ��ʾӦ����ֵ��
%  ���������
%     iStress --- Ӧ������ָʾ���������������ֵ
%                 1  --  x������Ӧ��
%                 2  --  y������Ӧ��
%                 3  --  ��Ӧ��
%                 4  --  ��һ��Ӧ��
%                 5  --  �ڶ���Ӧ��
%                 6  --  ����Ӧ��
%     nContour -- ��ֵ�ߵ�����
%     color  ---- ��ֵ����ɫ
%  ����ֵ��
%     ��
    global gNode gElement gNodeStress
   
    [element_number, dummy] = size( gElement ) ;
    [node_number, dummy] = size( gNode ) ;
    
    stressMin = min( gNodeStress( :, iStress ) ) ;
    stressMax = max( gNodeStress( :, iStress ) ) ;
    stressDelta = (stressMax-stressMin)/( nContour+1 ) ;
    stressContour = stressMin+stressDelta : stressDelta : stressMax - stressDelta ;
    for ie=1:1:element_number
        x = [ gNode( gElement( ie, 1 ), 1 ) ;
              gNode( gElement( ie, 2 ), 1 ) ;
              gNode( gElement( ie, 3 ), 1 ) ] ;
        y = [ gNode( gElement( ie, 1 ), 2 ) ;
              gNode( gElement( ie, 2 ), 2 ) ;
              gNode( gElement( ie, 3 ), 2 ) ] ;
        s = [ gNodeStress( gElement( ie, 1 ), iStress ) ;
              gNodeStress( gElement( ie, 2 ), iStress ) ;
              gNodeStress( gElement( ie, 3 ), iStress ) ] ;
        for is = 1:1:nContour  
            [smax, ismax] = max( s ) ;
            [smin, ismin] = min( s ) ;
            if stressContour(is) > smax | stressContour(is) < smin
                continue ;
            end
            
            x1 = x(ismin) + ( stressContour(is)- smin ) / (smax-smin) * ( x(ismax) - x(ismin) ) ;
            y1 = y(ismin) + ( stressContour(is)- smin ) / (smax-smin) * ( y(ismax) - y(ismin) ) ;
            
            for ismed=1:1:3
                if ismed ~= ismax & ismed ~= ismin
                    break ;
                end
            end
            
            if stressContour(is) < s( ismed )
                x2 = x(ismin) + (stressContour(is)-smin)/(s(ismed)-smin)*(x(ismed)-x(ismin)) ;
                y2 = y(ismin) + (stressContour(is)-smin)/(s(ismed)-smin)*(y(ismed)-y(ismin)) ;
            else
                x2 = x(ismed) + (stressContour(is)-s(ismed))/(smax-s(ismed))*(x(ismax)-x(ismed)) ;
                y2 = y(ismed) + (stressContour(is)-s(ismed))/(smax-s(ismed))*(y(ismax)-y(ismed)) ;
            end
          
            set( line( [x1;x2], [y1;y2] ), 'color', color ) ;
        end
    end
return

function PlotDisplacement( iDisp )
%  ��ʾλ����ͼ
%  ���������
%     iDisp --- λ�Ʒ���ָʾ���������������ֵ
%                 1  --  x����λ��
%                 2  --  y����λ��
%  ����ֵ��
%     ��
    global gNode gElement gDelta
    
    switch iDisp
    case 1
        title = ' x ����λ��' ;
    case 2
        title = ' y ����λ��' ;
    otherwise
        disp( sprintf( 'λ�Ʒ�������' ) ) ;
        return ;
    end
    
    figure ;
    axis equal ;
    axis off ;
    set( gcf, 'NumberTitle', 'off' ) ;
    set( gcf, 'Name', title ) ;
    
    dispMin = min( gDelta( iDisp:2:length(gDelta) ) ) ;
    dispMax = max( gDelta( iDisp:2:length(gDelta) ) ) ;
    caxis( [dispMin, dispMax] ) ;
    colormap( 'jet' ) ;

    [element_number, dummy] = size( gElement ) ;
    for ie=1:1:element_number
        x = [ gNode( gElement( ie, 1 ), 1 ) ;
              gNode( gElement( ie, 2 ), 1 ) ;
              gNode( gElement( ie, 3 ), 1 ) ] ;
        y = [ gNode( gElement( ie, 1 ), 2 ) ;
              gNode( gElement( ie, 2 ), 2 ) ;
              gNode( gElement( ie, 3 ), 2 ) ] ;
        c = [ gDelta( ( gElement( ie, 1 ) - 1 ) * 2 + iDisp ) ;
              gDelta( ( gElement( ie, 2 ) - 1 ) * 2 + iDisp ) ;
              gDelta( ( gElement( ie, 3 ) - 1 ) * 2 + iDisp ) ] ;
        set( patch( x, y, c ), 'EdgeColor', 'interp' ) ;
    end
    
    yTick = dispMin:(dispMax-dispMin)/10:dispMax ;
    Label = cell( 1, length(yTick) ); 
    for i=1:length(yTick)
        Label{i} = sprintf( '%.2e', yTick(i) ) ;
    end
    set( colorbar( 'vert' ), 'YTick', yTick, 'YTickLabelMode', 'Manual', 'YTickLabel', Label ) ;
    
    PlotDisplacementContour( iDisp, 10, 'white' ) ;
return

function PlotDisplacementContour( iDisp, nContour, color )
%  ��ʾλ�Ƶ�ֵ��
%  ���������
%     iDisp ----- λ�Ʒ���ָʾ���������������ֵ
%                 1  --  x����λ��
%                 2  --  y����λ��
%     nContour -- ��ֵ�ߵ�����
%     color  ---- ��ֵ����ɫ
%  ����ֵ��
%     ��
    global gNode gElement gDelta
   
    [element_number, dummy] = size( gElement ) ;
    [node_number, dummy] = size( gNode ) ;
    
    dispMin = min( gDelta( iDisp:2:length(gDelta) ) ) ;
    dispMax = max( gDelta( iDisp:2:length(gDelta) ) ) ;
    dispDelta = (dispMax-dispMin)/( nContour+1 ) ;
    dispContour = dispMin+dispDelta : dispDelta : dispMax - dispDelta ;
    for ie=1:1:element_number
        x = [ gNode( gElement( ie, 1 ), 1 ) ;
              gNode( gElement( ie, 2 ), 1 ) ;
              gNode( gElement( ie, 3 ), 1 ) ] ;
        y = [ gNode( gElement( ie, 1 ), 2 ) ;
              gNode( gElement( ie, 2 ), 2 ) ;
              gNode( gElement( ie, 3 ), 2 ) ] ;
        s = [ gDelta( ( gElement( ie, 1 ) - 1 ) * 2 + iDisp ) ;
              gDelta( ( gElement( ie, 2 ) - 1 ) * 2 + iDisp ) ;
              gDelta( ( gElement( ie, 3 ) - 1 ) * 2 + iDisp ) ] ;
        for is = 1:1:nContour  
            [smax, ismax] = max( s ) ;
            [smin, ismin] = min( s ) ;
            if dispContour(is) > smax | dispContour(is) < smin
                continue ;
            end
            
            x1 = x(ismin) + ( dispContour(is)- smin ) / (smax-smin) * ( x(ismax) - x(ismin) ) ;
            y1 = y(ismin) + ( dispContour(is)- smin ) / (smax-smin) * ( y(ismax) - y(ismin) ) ;
            
            for ismed=1:1:3
                if ismed ~= ismax & ismed ~= ismin
                    break ;
                end
            end
            
            if dispContour(is) < s( ismed )
                x2 = x(ismin) + (dispContour(is)-smin)/(s(ismed)-smin)*(x(ismed)-x(ismin)) ;
                y2 = y(ismin) + (dispContour(is)-smin)/(s(ismed)-smin)*(y(ismed)-y(ismin)) ;
            else
                x2 = x(ismed) + (dispContour(is)-s(ismed))/(smax-s(ismed))*(x(ismax)-x(ismed)) ;
                y2 = y(ismed) + (dispContour(is)-s(ismed))/(smax-s(ismed))*(y(ismax)-y(ismed)) ;
            end
          
            set( line( [x1;x2], [y1;y2] ), 'color', color ) ;
        end
    end
return