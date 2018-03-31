function exam4_2( file_in )
% ������Ϊ�����µĵڶ������������������ε�Ԫ������������������µı��κ�Ӧ��
%      exam4_2( filename )
%  ��������� 
%      file_in  ---------- ����Ԫģ���ļ�

% ����ȫ�ֱ���
%      gNode ------------- �ڵ�����
%      gElement ---------- ��Ԫ����
%      gMaterial --------- ��������
%      gBC1 -------------- ��һ��Լ������
%      gK ---------------- ����նȾ���
%      gDelta ------------ ����ڵ�����
%      gNodeStress ------- �ڵ�Ӧ��
%      gElementStress ---- ��ԪӦ��
    global gNode gElement gMaterial gBC1 gK gDelta gNodeStress gElementStress

    if nargin < 1
        file_in = 'exam4_2.dat' ;
    end
    
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

    % ��������Ԫģ�Ͳ���⣬������
    FemModel( file_in ) ;          % ��������Ԫģ��
    SolveModel ;                   % �������Ԫģ��
    SaveResults( file_out ) ;      % ���������
    
    % �������
    disp( sprintf( '������������������������ļ� %s ��', file_out ) ) ;
    disp( sprintf( '����ʹ�ú������ exam4_2_post.m ��ʾ������' ) ) ;
return ;

function FemModel(filename)
%  ��������Ԫģ��
%  ���������
%      filename --- ����Ԫģ���ļ�
%  ����ֵ��
%      ��
%  ˵����
%      �ú�������ƽ�����������Ԫģ�����ݣ�
%        gNode ------- �ڵ㶨��
%        gElement ---- ��Ԫ����
%        gMaterial --- ���϶��壬��������ģ�������Ľ���������Ŀ�����Ծ�
%        gBC1 -------- Լ������

    global gNode gElement gMaterial gBC1
    
    % ���ļ�
    fid = fopen( filename, 'r' ) ;
    
    % ��ȡ�ڵ�����
    node_number = fscanf( fid, '%d', 1 ) ;
    gNode = zeros( node_number, 2 ) ;
    for i=1:node_number
        dummy = fscanf( fid, '%d', 1 ) ;
        gNode( i, : ) = fscanf( fid, '%f', [1, 2] ) ;
    end
    
    % ��ȡ��Ԫ����
    element_number = fscanf( fid, '%d', 1 ) ;
    gElement = zeros( element_number, 4 ) ;
    for i=1:element_number
        dummy = fscanf( fid, '%d', 1 ) ;
        gElement( i, : ) = fscanf( fid, '%d', [1, 4] ) ;
    end
    
    % ��ȡ������Ϣ
    material_number = fscanf( fid, '%d', 1 ) ;
    gMaterial = zeros( material_number, 4 ) ;
    for i=1:material_number
        dummy = fscanf( fid, '%d', 1 ) ;
        gMaterial( i, : ) = fscanf( fid, '%f', [1,4] ) ;
    end
    
    % ��ȡ�߽�����
    bc1_number = fscanf( fid, '%d', 1 ) ;
    gBC1 = zeros( bc1_number, 3 ) ;
    for i=1:bc1_number
        gBC1( i, 1 ) = fscanf( fid, '%d', 1 ) ;
        gBC1( i, 2 ) = fscanf( fid, '%d', 1 ) ;
        gBC1( i, 3 ) = fscanf( fid, '%f', 1 ) ;
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
%        5. ���㵥ԪӦ���ͽڵ�Ӧ��

    global gNode gElement gMaterial gBC1 gK gDelta gNodeStress gElementStress

    % step1. ��������նȾ���ͽڵ�������
    [node_number,dummy] = size( gNode ) ;
    gK = sparse( node_number * 2, node_number * 2 ) ;
    f = sparse( node_number * 2, 1 ) ;

    % step2. ���㵥Ԫ�նȾ��󣬲����ɵ�����նȾ�����
    [element_number,dummy] = size( gElement ) ;
    for ie=1:1:element_number
        disp( sprintf(  '����նȾ��󣬵�ǰ��Ԫ: %d', ie  ) ) ;
        k = StiffnessMatrix( ie ) ;
        AssembleStiffnessMatrix( ie, k ) ;
    end

    % step3. �������ز����ĵ�Ч�ڵ���
    for ie=1:1:element_number
        disp( sprintf(  '�������صĵ�Ч�ڵ�������ǰ��Ԫ: %d', ie  ) ) ;
        egf = EquivalentGravityForce( ie ) ;
        i = gElement( ie, 1 ) ;
        j = gElement( ie, 2 ) ;
        m = gElement( ie, 3 ) ;
        f( (i-1)*2+1 : (i-1)*2+2 ) = f( (i-1)*2+1 : (i-1)*2+2 ) + egf( 1:2 ) ;
        f( (j-1)*2+1 : (j-1)*2+2 ) = f( (j-1)*2+1 : (j-1)*2+2 ) + egf( 3:4 ) ;
        f( (m-1)*2+1 : (m-1)*2+2 ) = f( (m-1)*2+1 : (m-1)*2+2 ) + egf( 5:6 ) ;
    end
  
    % step4. ����Լ���������޸ĸնȾ���ͽڵ������������ó˴�����
    [bc_number,dummy] = size( gBC1 ) ;
    for ibc=1:1:bc_number
        n = gBC1(ibc, 1 ) ;
        d = gBC1(ibc, 2 ) ;
        m = (n-1)*2 + d ;
        f(m) = gBC1(ibc, 3)* gK(m,m) * 1e15 ;
        gK(m,m) = gK(m,m) * 1e15 ;
    end

    % step 5. ��ⷽ���飬�õ��ڵ�λ������
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
        S = zeros( 1, 3 ) ;                         
        A = 0 ;
        for ie=1:1:element_number
            for k=1:1:3
                if i == gElement( ie, k ) 
                    area= ElementArea( ie ) ;
                    S = S + gElementStress(ie,1:3 ) * area ;
                    A = A + area ;
                    break ;
                end
            end
        end
        gNodeStress(i,1:3) = S / A ;
        gNodeStress(i,6) = 0.5*sqrt( (gNodeStress(i,1)-gNodeStress(i,2))^2 + 4*gNodeStress(i,3)^2 ) ;
        gNodeStress(i,4) = 0.5*(gNodeStress(i,1)+gNodeStress(i,2)) + gNodeStress(i,6) ;
        gNodeStress(i,5) = 0.5*(gNodeStress(i,1)+gNodeStress(i,2)) - gNodeStress(i,6) ;
    end
return

function k = StiffnessMatrix( ie )
%  ���㵥Ԫ�նȾ���
%  �������:
%     ie ----  ��Ԫ��
%  ����ֵ:
%     k  ----  ��Ԫ�նȾ���

    global gNode gElement gMaterial 
    k = zeros( 6, 6 ) ;
    E  = gMaterial( gElement(ie, 4), 1 ) ;
    mu = gMaterial( gElement(ie, 4), 2 ) ;
    h  = gMaterial( gElement(ie, 4), 3 ) ;
    xi = gNode( gElement( ie, 1 ), 1 ) ;
    yi = gNode( gElement( ie, 1 ), 2 ) ;
    xj = gNode( gElement( ie, 2 ), 1 ) ;
    yj = gNode( gElement( ie, 2 ), 2 ) ;
    xm = gNode( gElement( ie, 3 ), 1 ) ;
    ym = gNode( gElement( ie, 3 ), 2 ) ;
    ai = xj*ym - xm*yj ;
    aj = xm*yi - xi*ym ;
    am = xi*yj - xj*yi ;
    bi = yj - ym ;
    bj = ym - yi ;
    bm = yi - yj ;
    ci = -(xj-xm) ;
    cj = -(xm-xi) ;
    cm = -(xi-xj) ;
    area = (ai+aj+am)/2 ;
    B = [bi  0 bj  0 bm  0
          0 ci  0 cj  0 cm
         ci bi cj bj cm bm] ;
    B = B/2/area ;
    D = [ 1-mu    mu      0
           mu    1-mu     0
            0      0   (1-2*mu)/2] ;
    D = D*E/(1-2*mu)/(1+mu) ;
    k = transpose(B)*D*B*h*abs(area) ;    
return

function B = MatrixB( ie )
%  ���㵥Ԫ��Ӧ�����B
%  �������:
%     ie ----  ��Ԫ��
%  ����ֵ:
%     B  ----  ��ԪӦ�����
    global gNode gElement
    xi = gNode( gElement( ie, 1 ), 1 ) ;
    yi = gNode( gElement( ie, 1 ), 2 ) ;
    xj = gNode( gElement( ie, 2 ), 1 ) ;
    yj = gNode( gElement( ie, 2 ), 2 ) ;
    xm = gNode( gElement( ie, 3 ), 1 ) ;
    ym = gNode( gElement( ie, 3 ), 2 ) ;
    ai = xj*ym - xm*yj ;
    aj = xm*yi - xi*ym ;
    am = xi*yj - xj*yi ;
    bi = yj - ym ;
    bj = ym - yi ;
    bm = yi - yj ;
    ci = -(xj-xm) ;
    cj = -(xm-xi) ;
    cm = -(xi-xj) ;
    area = (ai+aj+am)/2 ;
    B = [bi  0 bj  0 bm  0
          0 ci  0 cj  0 cm
         ci bi cj bj cm bm] ;
    B = B/2/area ;
return

function area = ElementArea( ie )
%  ���㵥Ԫ���
%  �������:
%     ie ----  ��Ԫ��
%  ����ֵ:
%     area  ----  ��Ԫ���
    global gNode gElement gMaterial 

    xi = gNode( gElement( ie, 1 ), 1 ) ;
    yi = gNode( gElement( ie, 1 ), 2 ) ;
    xj = gNode( gElement( ie, 2 ), 1 ) ;
    yj = gNode( gElement( ie, 2 ), 2 ) ;
    xm = gNode( gElement( ie, 3 ), 1 ) ;
    ym = gNode( gElement( ie, 3 ), 2 ) ;
    ai = xj*ym - xm*yj ;
    aj = xm*yi - xi*ym ;
    am = xi*yj - xj*yi ;
    area = abs((ai+aj+am)/2) ;
return

function AssembleStiffnessMatrix( ie, k )
%  �ѵ�Ԫ�նȾ��󼯳ɵ�����նȾ���
%  �������:
%      ie  --- ��Ԫ��
%      k   --- ��Ԫ�նȾ���
%  ����ֵ:
%      ��
    global gElement gK
    for i=1:1:3
        for j=1:1:3
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

function egf = EquivalentGravityForce( ie )
%  ���㵥Ԫ���صĵ�Ч�ڵ���
%  �������
%      ie  ----- ��Ԫ��
%  ����ֵ
%      egf ----- ���صĵ�Ч�ڵ��� 
    global gElement gMaterial
    
    area = ElementArea( ie ) ;
    h    = gMaterial( gElement( ie, 4 ), 3 ) ;
    ro   = gMaterial( gElement( ie, 4 ), 4 ) ;
    w    = area * h * ro ;
    egf  = -w/3 * [0; 1; 0; 1; 0; 1] ;
return

function es = ElementStress( ie )
%  ���㵥Ԫ��Ӧ������
%  �������
%      ie  ----- ��Ԫ��
%  ����ֵ
%      es ----- ��ԪӦ����������1��6���� [sx, sy, txy, s1, s2, tmax]
    global gElement gDelta  gMaterial
    
    es = zeros( 1, 6 ) ;   % ��Ԫ��Ӧ������
    de = zeros( 6, 1 ) ;   % ��Ԫ�ڵ�λ������
    E  = gMaterial( gElement(ie, 4), 1 ) ;
    mu = gMaterial( gElement(ie, 4), 2 ) ;
    D = [ 1-mu    mu      0
           mu    1-mu     0
            0      0   (1-2*mu)/2] ;
    D = D*E/(1-2*mu)/(1+mu) ;
    B = MatrixB( ie ) ;
    for j=1:1:3
        de( 2*j-1 ) = gDelta( 2*gElement( ie, j )-1 ) ;
        de( 2*j   ) = gDelta( 2*gElement( ie, j )   ) ;
    end
    es(1:3) = D * B * de ;
    es(6) = 0.5*sqrt((es(1)-es(2))^2 + 4*es(3)^2 ) ;
    es(4) = 0.5*(es(1)+es(2)) + es(6) ;
    es(5) = 0.5*(es(1)+es(2)) - es(6) ;
return

function SaveResults( file_out )
%  ��ʾ������
%  ���������
%     file_out  --- �������ļ�
%  ����ֵ��
%     ��

    global gNode gElement gMaterial gBC1 gDelta gNodeStress gElementStress
    save( file_out, 'gNode', 'gElement', 'gMaterial', 'gBC1', 'gDelta', 'gNodeStress', 'gElementStress' ) ;
return