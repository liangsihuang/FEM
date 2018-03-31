function exam4_2_post( file_in )
%  三角形单元的后处理程序
%  输入参数
%     file_in  -----  计算结果文件，由 exam4_2.m 生成
%  返回值
%     无
    global gNode gElement gMaterial gBC1 gDelta gNodeStress gElementStress
    
    if nargin < 1 
        file_in = 'exam4_2.mat' ;
    end
    
    load( file_in ) ;
    gDelta = full( gDelta ) ;
    DisplayModel ;            %  显示初始模型
    for i=1:6
        PlotStress( i ) ;     %  显示应力云图    
    end
    PlotDisplacement( 1 ) ;   %  显示x方向位移云图
    PlotDisplacement( 2 ) ;   %  显示y方向位移云图
return

function DisplayModel
%  用图形方式显示有限元模型
%  输入参数：
%     无
%  返回值：
%     无
    global gNode gElement gMaterial gBC1 
    
    figure ;
    axis equal ;
    axis off ;
    set( gcf, 'NumberTitle', 'off' ) ;
    set( gcf, 'Name', '有限元模型' ) ;
    
    % 根据不同的材料，显示单元颜色
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
    
    % 显示边界条件
    DisplayBC( 'blue' ) ;
    
return

function DisplayBC( color )
%  用图形方式显示有限元模型的边界条件
%  输入参数：
%     color  ----  边界条件的颜色
%  返回值：
%     无
    global gNode gBC1
    
    % 确定边界条件的大小
    xmin = min( gNode(:,1) ) ;
    xmax = max( gNode(:,1) ) ;
    factor = ( xmax - xmin ) / 25 ;
    
    [bc1_number,dummy] = size( gBC1 ) ;
    dBCSize = factor ;
    for i=1:bc1_number
        if( gBC1( i, 2 ) == 1 )  %  x方向约束
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
        else                      %  y方向约束
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
%  显示应力云图
%  输入参数：
%     iStress --- 应力分量指示，它可以是下面的值
%                 1  --  x方向正应力
%                 2  --  y方向正应力
%                 3  --  剪应力
%                 4  --  第一主应力
%                 5  --  第二主应力
%                 6  --  最大剪应力
%  返回值：
%     无
    global gNode gElement gNodeStress
    
    switch iStress
    case 1
        title = ' x 方向正应力' ;
    case 2
        title = ' y 方向正应力' ;
    case 3
        title = ' 剪应力' ;
    case 4
        title = ' 最大主应力' ;
    case 5
        title = ' 最小主应力' ;
    case 6
        title = ' 最大剪应力' ;
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
%  显示应力等值线
%  输入参数：
%     iStress --- 应力分量指示，它可以是下面的值
%                 1  --  x方向正应力
%                 2  --  y方向正应力
%                 3  --  剪应力
%                 4  --  第一主应力
%                 5  --  第二主应力
%                 6  --  最大剪应力
%     nContour -- 等值线的条数
%     color  ---- 等值线颜色
%  返回值：
%     无
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
%  显示位移云图
%  输入参数：
%     iDisp --- 位移分量指示，它可以是下面的值
%                 1  --  x方向位移
%                 2  --  y方向位移
%  返回值：
%     无
    global gNode gElement gDelta
    
    switch iDisp
    case 1
        title = ' x 方向位移' ;
    case 2
        title = ' y 方向位移' ;
    otherwise
        disp( sprintf( '位移分量错误' ) ) ;
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
%  显示位移等值线
%  输入参数：
%     iDisp ----- 位移分量指示，它可以是下面的值
%                 1  --  x方向位移
%                 2  --  y方向位移
%     nContour -- 等值线的条数
%     color  ---- 等值线颜色
%  返回值：
%     无
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