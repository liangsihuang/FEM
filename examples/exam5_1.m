function exam5_1
% ������Ϊ�����µĵ�һ�����������������嵥Ԫ���㵥������˵ĺ����Ӧ��
%  ��������� 
%      ��

% ����ȫ�ֱ���
%      gNode ------ �ڵ�����
%      gElement --- ��Ԫ����
%      gMaterial -- ��������
%      gBC1 ------- ��һ��Լ������
%      gDF -------- �ֲ���
%      gNF -------- ������
%      gK --------- ����նȾ���
%      gDelta ----- ����ڵ�����

    file_in = 'exam5_1.dat' ;
    
    % ����ļ��Ƿ����
    if exist( file_in ) == 0
        disp( sprintf( '�����ļ� %s ������', file_in ) )
        disp( sprintf( '������ֹ' ) )
        return ;
    end
    
    % ���������ļ�������������ļ�����
    [path_str,name_str,ext_str] = fileparts( file_in ) ;
    ext_str_out = '.mat' ;
    file_out = fullfile( path_str, [name_str, ext_str_out] ) ;
    
    % �������ļ��Ƿ����
    if exist( file_out ) ~= 0 
        answer = input( sprintf( '�ļ� %s �Ѿ����ڣ��Ƿ񸲸�? ( Yes / [No] ):  ', file_out ), 's' ) ;
        if length( answer ) == 0
            answer = 'no' ;
        end
        
        answer = lower( answer ) ;
        if answer ~= 'y' | answer ~= 'yes' 
            disp( sprintf( '��ʹ��������ļ������򱸷����е��ļ�' ) ) ;
            disp( sprintf( '������ֹ' ) );
            return ; 
        end
    end


    FemModel( file_in ) ;       % ��������Ԫģ��
    SolveModel ;                % �������Ԫģ��
    SaveResults( file_out ) ;   % ���������
    DisplayResults(file_out) ;  % ��ʾ���
return ;

function FemModel( file_in )
%  ��������Ԫģ��
%  ���������
%      file_in --- �����ļ�����
%  ����ֵ��
%  ˵����
    global gNode gElement gMaterial gBC1 gDF gNF
    
    % ���ļ�
    fid = fopen( file_in, 'r' ) ;
    
    % �ڵ�����
    node_number = fscanf( fid, '%d', 1 ) ;
    gNode = zeros( node_number, 3 ) ;
    for i = 1:node_number
        dummy = fscanf( fid, '%d', 1 ) ;
        gNode( i, : ) = fscanf( fid, '%f', [1 3] ) ;
    end
    
    % ��Ԫ����
    element_number = fscanf( fid, '%d', 1 ) ;
    gElement = zeros( element_number, 5 ) ;
    for i = 1:element_number
        dummy = fscanf( fid, '%d', 1 ) ;
        gElement( i, : ) = fscanf( fid, '%d', [1 5] ) ;
    end
    
    % ��������
    material_number = fscanf( fid, '%d', 1 ) ;
    gMaterial = zeros( material_number, 2 ) ;
    for i = 1:material_number
        dummy = fscanf( fid, '%d', 1 ) ;
        gMaterial( i, : ) = fscanf( fid, '%f', [1 2] ) ;
    end
    
    % ��һ��Լ������
    bc1_number = fscanf( fid, '%d', 1 ) ;
    gBC1 = zeros( bc1_number, 3 ) ;
    for i=1:bc1_number
        gBC1(i,1) = fscanf( fid, '%d', 1 ) ;  % ����
        gBC1(i,2) = fscanf( fid, '%d', 1 ) ;  % ���ɶȺ�
        gBC1(i,3) = fscanf( fid, '%f', 1 ) ;  % Լ��ֵ
    end
    
    % �ֲ�ѹ�������Էֲ���
    df_number = fscanf( fid, '%d', 1 ) ;
    gDF = zeros( df_number, 7 ) ;
    for i=1:df_number
        gDF(i,1:3) = fscanf( fid, '%d', [1 3] ) ;  % ��������
        gDF(i,4:6) = fscanf( fid, '%f', [1 3] ) ;  % ��������ϵ�ѹ��ֵ 
    end
    
    % ������
    nf_number = fscanf( fid, '%d', 1 ) ;
    gNF = zeros( nf_number, 3 ) ;
    for i=1:nf_number
        gNF(i,1) = fscanf( fid, '%d', 1 ) ;  % ����
        gNF(i,2) = fscanf( fid, '%d', 1 ) ;  % ���ɶȺ�
        gNF(i,3) = fscanf( fid, '%f', 1 ) ;  % ��������С
    end
    
    % �ر��ļ�
    fclose( fid ) ;
return

function SolveModel
%  �������Ԫģ��
%  ���������
%     ��
%  ����ֵ��
%     ��
%  ˵����
%      �ú����������Ԫģ�ͣ���������
%        1. ���㵥Ԫ�նȾ��󣬼�������նȾ���
%        2. ���㵥Ԫ�ĵ�Ч�ڵ�������������ڵ�������
%        3. ����Լ���������޸�����նȾ���ͽڵ�������
%        4. ��ⷽ���飬�õ�����ڵ�λ������

    global gNode gElement gMaterial gBC1 gDF gNF gK gDelta gNodeStress gElementStress

    % step1. ��������նȾ���ͽڵ�������
    [node_number,dummy] = size( gNode ) ;
    gK = sparse( node_number * 3, node_number * 3 ) ;
    f = sparse( node_number * 3, 1 ) ;

    % step2. ���㵥Ԫ�նȾ��󣬲����ɵ�����նȾ�����
    [element_number,dummy] = size( gElement ) ;
    for ie=1:1:element_number
        disp( sprintf(  '����ڵ�նȾ��󲢼��ɣ���ǰ��Ԫ: %d', ie  ) ) ;
        k = StiffnessMatrix( ie ) ;
        AssembleStiffnessMatrix( ie, k ) ;
    end

    % step3. ����ֲ�ѹ���ĵ�Ч�ڵ����������ɵ�����ڵ���������
    [df_number, dummy] = size( gDF ) ;
    for idf = 1:1:df_number
        enf = EquivalentDistPressure( gDF(idf,1:3), gDF(idf,4:6) ) ;
        i = gDF(idf, 1) ;
        j = gDF(idf, 2) ;
        m = gDF(idf, 3) ;
        f( (i-1)*3+1:(i-1)*3+3 ) = f( (i-1)*3+1:(i-1)*3+3 ) + enf( 1:3 ) ;
        f( (j-1)*3+1:(j-1)*3+3 ) = f( (j-1)*3+1:(j-1)*3+3 ) + enf( 4:6 ) ;
        f( (m-1)*3+1:(m-1)*3+3 ) = f( (m-1)*3+1:(m-1)*3+3 ) + enf( 7:9 ) ;
    end
  
    % step4. �Ѽ��������ɵ�����ڵ���������
    [nf_number, dummy] = size( gNF ); 
    for i = 1:1:nf_number
        in = gNF( i, 1 ) ;
        id = gNF( i, 2 ) ;
        f( (in-1)*3 + id ) = f( (in-1)*3 + id ) + gNF( i, 3 ) ;
    end
    
    % step5. ����Լ���������޸ĸնȾ���ͽڵ������������ó˴�����
    [bc_number,dummy] = size( gBC1 ) ;
    for ibc=1:1:bc_number
        n = gBC1(ibc, 1 ) ;
        d = gBC1(ibc, 2 ) ;
        m = (n-1)*3 + d ;
        f(m) = gBC1(ibc, 3)* gK(m,m) * 1e15 ;
        gK(m,m) = gK(m,m) * 1e15 ;
    end

    % step 6. ��ⷽ���飬�õ��ڵ�λ������
    gDelta = gK \ f ;
    
    % step 6. ���㵥ԪӦ��
    gElementStress = zeros( element_number, 6 ) ;
    for ie=1:element_number 
        disp( sprintf(  '���㵥ԪӦ������ǰ��Ԫ: %d', ie  ) ) ;
        es = ElementStress( ie ) ;
        gElementStress( ie, : ) = es ;
    end
    
    % step 7. ����ڵ�Ӧ��(�����ƽڵ��Ȩƽ��)
    gNodeStress = zeros( node_number, 6 ) ;       
    for i=1:node_number                         
        disp( sprintf(  '����ڵ�Ӧ������ǰ���: %d', i  ) ) ;
        S = zeros( 1, 6 ) ;                         
        V = 0 ;
        for ie=1:1:element_number
            for k=1:1:4
                if i == gElement( ie, k ) 
                    vol = ElementVolume( ie ) ;
                    S = S + gElementStress(ie,:) * vol ;
                    V = V + vol ;
                    break ;
                end
            end
        end
        gNodeStress(i,:) = S / V ;
    end
return

function k = StiffnessMatrix( ie )
%  ���㵥Ԫ�նȾ���
%  �������:
%     ie ----  ��Ԫ��
%  ����ֵ:
%     k  ----  ��Ԫ�նȾ���

    global gNode gElement gMaterial 
    
    % ��ȡ�������
    x = gNode( gElement( ie, : ), 1 ) ;
    y = gNode( gElement( ie, : ), 2 ) ;
    z = gNode( gElement( ie, : ), 3 ) ;
    
    % ���� 6V
    V6 = det( [1 x(1) y(1) z(1)
               1 x(2) y(2) z(2)
               1 x(3) y(3) z(3)
               1 x(4) y(4) z(4) ] ) ;
           
    if V6 < 0
        disp( sprintf( '���棺��Ԫ %d �Ľ������˳��������', ie ) ) ;
        pause ;
    end
    % ����Ӧ�����
    B = zeros( 6, 12 ) ;
    for i=1:4
        j = mod(i,  4) + 1 ;
        m = mod(i+1,4) + 1 ;
        p = mod(i+2,4) + 1 ;
        bi = - det( [ 1 y(j) z(j)
                      1 y(m) z(m)
                      1 y(p) z(p) ] ) ;
        ci =   det( [ 1 x(j) z(j)
                      1 x(m) z(m)
                      1 x(p) z(p) ] ) ;
        di = - det( [ 1 x(j) y(j)
                      1 x(m) y(m)
                      1 x(p) y(p) ] ) ;
        B(:,(i-1)*3+1:(i-1)*3+3 ) = (-1)^(i+1)*[bi   0   0
                                                 0  ci   0
                                                 0   0  di
                                                 0  di  ci
                                                di   0  bi
                                                ci  bi   0] ;
    end
    B = B / V6 ;
    
    % ���㵯�Ծ���
    E = gMaterial( gElement( ie, 5 ), 1 ) ;
    mu = gMaterial( gElement( ie, 5 ), 2 ) ;
    A1 = mu/(1-mu) ;
    A2 = (1-2*mu)/2/(1-mu) ;
    A3 = E*(1-mu)/(1+mu)/(1-2*mu) ;
    D = [  1  A1  A1   0   0   0
          A1   1  A1   0   0   0
          A1  A1   1   0   0   0
           0   0   0  A2   0   0 
           0   0   0  0   A2   0
           0   0   0  0    0  A2] ;
    D = D * A3 ;
    
    % ���㵥Ԫ�նȾ���
    V = abs(V6)/6 ;
    k = transpose( B ) * D * B * V ;
return

function AssembleStiffnessMatrix( ie, k )
%  �ѵ�Ԫ�նȾ��󼯳ɵ�����նȾ���
%  �������:
%      ie  --- ��Ԫ��
%      k   --- ��Ԫ�նȾ���
%  ����ֵ:
%      ��
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

function enf = EquivalentDistPressure( node, pressure )
%   �������Էֲ�ѹ���ĵ�Ч�ڵ���
%   �������:
%      node ----------  ����
%      pressure ------  �����Ŷ�Ӧ��ѹ��ֵ
%   ����ֵ:
%      enf -----------  ��Ч�ڵ�������
    global gNode
    
    enf = zeros( 9, 1 ) ;
    % ��������ѹ���������ε����
    x = gNode( node, 1 ) ;
    y = gNode( node, 2 ) ;
    z = gNode( node, 3 ) ;
    b = - det( [1  y(1) z(1)
                1  y(2) z(2)
                1  y(3) z(3)] ) ;
    c =   det( [1  x(1) z(1)
                1  x(2) z(2)
                1  x(3) z(3)] ) ;
    d = - det( [1  x(1) y(1)
                1  x(2) y(2)
                1  x(3) y(3)] ) ;
    A2 = sqrt( b^2 + c^2 + d^2 ) ;
    A = A2 / 2 ;

    % ������������ϵĵ�Ч�ڵ���
    f1 = ( pressure(1)   + pressure(2)/2 + pressure(3)/2 ) / 6 * A ;
    f2 = ( pressure(1)/2 + pressure(2)   + pressure(3)/2 ) / 6 * A ;
    f3 = ( pressure(1)/2 + pressure(2)/2 + pressure(3)   ) / 6 * A ;
    
    % �������������ϵķ���
    nx = -b/A2 ;
    ny = -c/A2 ;
    nz = -d/A2 ;
    enf(1:3) = f1*[nx;ny;nz] ;
    enf(4:6) = f2*[nx;ny;nz] ;
    enf(7:9) = f3*[nx;ny;nz] ;
return

function es = ElementStress( ie )
%  ���㵥Ԫ��Ӧ������
%  �������
%      ie  ----- ��Ԫ��
%  ����ֵ
%      es ----- ��ԪӦ����������1��6���� [sx, sy, sz, tyz, txz, txy]
    global gElement gDelta gMaterial
    
    de = zeros( 12, 1 ) ;   % ��Ԫ�ڵ�λ������
    E  = gMaterial( gElement(ie, 5), 1 ) ;
    mu = gMaterial( gElement(ie, 5), 2 ) ;
    A1 = mu/(1-mu) ;
    A2 = (1-2*mu)/2/(1-mu) ;
    A3 = E*(1-mu)/(1+mu)/(1-2*mu) ;
    D = [  1  A1  A1   0   0   0
          A1   1  A1   0   0   0
          A1  A1   1   0   0   0
           0   0   0  A2   0   0 
           0   0   0   0  A2   0
           0   0   0   0   0  A2] ;
    D = D * A3 ;
    B = MatrixB( ie ) ;
    for j=1:1:4
        de( 3*j-2 ) = gDelta( 3*gElement( ie, j )-2 ) ;
        de( 3*j-1 ) = gDelta( 3*gElement( ie, j )-1 ) ;
        de( 3*j   ) = gDelta( 3*gElement( ie, j )   ) ;
    end
    es = D * B * de ;
    es = transpose( es ) ;
return

function vol = ElementVolume( ie )
%  ���㵥Ԫ�����
%  �������
%      ie  ----- ��Ԫ��
%  ����ֵ
%      vol ----- ��Ԫ���

    global gNode gElement
    
    % ��ȡ�������
    x = gNode( gElement( ie, : ), 1 ) ;
    y = gNode( gElement( ie, : ), 2 ) ;
    z = gNode( gElement( ie, : ), 3 ) ;
    
    % ���� 6V
    V6 = det( [1 x(1) y(1) z(1)
               1 x(2) y(2) z(2)
               1 x(3) y(3) z(3)
               1 x(4) y(4) z(4) ] ) ;
    
    vol = abs( V6/6 ) ;    
return

function B = MatrixB(ie)
%  ���㵥Ԫ��Ӧ�����B
%  �������:
%     ie ----  ��Ԫ��
%  ����ֵ:
%     B  ----  ��ԪӦ�����

    global gNode gElement 
    
    % ��ȡ�������
    x = gNode( gElement( ie, : ), 1 ) ;
    y = gNode( gElement( ie, : ), 2 ) ;
    z = gNode( gElement( ie, : ), 3 ) ;
    
    % ���� 6V
    V6 = det( [1 x(1) y(1) z(1)
               1 x(2) y(2) z(2)
               1 x(3) y(3) z(3)
               1 x(4) y(4) z(4) ] ) ;
           
    % ����Ӧ�����
    B = zeros( 6, 12 ) ;
    for i=1:4
        j = mod(i,  4) + 1 ;
        m = mod(i+1,4) + 1 ;
        p = mod(i+2,4) + 1 ;
        bi = - det( [ 1 y(j) z(j)
                      1 y(m) z(m)
                      1 y(p) z(p) ] ) ;
        ci =   det( [ 1 x(j) z(j)
                      1 x(m) z(m)
                      1 x(p) z(p) ] ) ;
        di = - det( [ 1 x(j) y(j)
                      1 x(m) y(m)
                      1 x(p) y(p) ] ) ;
        B(:,(i-1)*3+1:(i-1)*3+3 ) = (-1)^(i+1)*[bi   0   0
                                                 0  ci   0
                                                 0   0  di
                                                 0  di  ci
                                                di   0  bi
                                                ci  bi   0] ;
    end
    B = B / V6 ;
return

function SaveResults( file_out ) 
%  ���������
%  ���������
%     file_out ---- �����ļ���
%  ����ֵ��
%     ��

    global gNode gElement gMaterial gBC1 gDF gNF gK gDelta gNodeStress gElementStress
    save( file_out, 'gNode', 'gElement', 'gMaterial', 'gBC1', ...
          'gDF', 'gNF', 'gDelta', 'gNodeStress', 'gElementStress' ) ;
return

function DisplayResults(file_out)
%  ��ʾ������
%  ���������
%     ��
%  ����ֵ��
%     ��
    y = [3,3.5,4,4.5] ;
    color = 'bgrc' ;
    [x,sy,syi] = exam5_1_post( file_out,y,0:0.25:4.5 ) ;
    figure ;
    hold ;
    for i=1:length(y)
        plot( x, sy(i,:),color(i) ) ;
    end
    legend( 'y=3.0','y=3.5','y=4.0','y=4.5' ) ;
    title( '�뼯������ͬ���봦����Ӧ���غ����ķֲ�' ) ;
    xlabel( 'x' ) ;
    ylabel( '��Ӧ��' ) ;
    hold off ;
return