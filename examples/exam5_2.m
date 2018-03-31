function exam5_2
% ������Ϊ�����µĵڶ���������������������ԳƵ�Ԫ����Բ�δ�ֱ���������µػ��ı���

% ����ȫ�ֱ���
%      gNode ------------- �ڵ�����
%      gElement ---------- ��Ԫ����
%      gMaterial --------- ��������
%      gBC1 -------------- ��һ��Լ������
%      gK ---------------- ����նȾ���
%      gDelta ------------ ����ڵ�����
    global gNode gElement gMaterial gBC1 gK gDelta

    file_in = 'exam5_2.dat' ;
    
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
    DisplayResults ;               % ��ʾ������    
return ;

function FemModel(filename)
%  ��������Ԫģ��
%  ���������
%      filename --- ����Ԫģ���ļ�
%  ����ֵ��
%      ��
%  ˵����
%      �ú�����������Ԫģ�����ݣ�
%        gNode ------- �ڵ㶨��
%        gElement ---- ��Ԫ����
%        gMaterial --- ���϶��壬��������ģ�������Ľ���������Ŀ�����Ծ�
%        gBC1--------- ��һ��Լ������
%        gDFStep ----- �ֲ������ز�
%        gDF --------- �ֲ���

    global gNode gElement gMaterial gBC1 gDF gDFStep
    
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
    gMaterial = zeros( material_number, 2 ) ;
    for i=1:material_number
        dummy = fscanf( fid, '%d', 1 ) ;
        gMaterial( i, : ) = fscanf( fid, '%f', [1,2] ) ;
    end
    
    % ��ȡ�߽�����
    bc1_number = fscanf( fid, '%d', 1 ) ;
    gBC1 = zeros( bc1_number, 3 ) ;
    for i=1:bc1_number
        gBC1( i, 1 ) = fscanf( fid, '%d', 1 ) ;
        gBC1( i, 2 ) = fscanf( fid, '%d', 1 ) ;
        gBC1( i, 3 ) = fscanf( fid, '%f', 1 ) ;
    end
    
    % ��ȡ�ֲ�ѹ������
    dfstep_number = fscanf( fid, '%d', 1 ) ;
    gDFStep = fscanf( fid, '%d', [dfstep_number,1] ) ;
    
    % ��ȡ�ֲ�ѹ��
    gDF = [] ;
    for idfstep = 1:dfstep_number
        df_number = fscanf( fid, '%d', 1 ) ;
        DF = zeros( df_number, 4 ) ;
        for i=1:df_number
            DF( i, 1:2 ) = fscanf( fid, '%d', [1, 2] ) ;   % ����
            DF( i, 3:4 ) = fscanf( fid, '%f', [1, 2] ) ;   % ��Ӧ�ڽ��ŵ�ѹ��ֵ
        end
        gDF = [ gDF; DF ] ;
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

    global gNode gElement gMaterial gBC1 gDF gDFStep gK gDelta

    % step1. ��������նȾ���ͽڵ�������
    [node_number,dummy] = size( gNode ) ;
    gK = sparse( node_number * 2, node_number * 2 ) ;
    f = sparse( node_number * 2, length(gDFStep) ) ;

    % step2. ���㵥Ԫ�նȾ��󣬲����ɵ�����նȾ�����
    [element_number,dummy] = size( gElement ) ;
    for ie=1:1:element_number
        disp( sprintf(  '����նȾ��󣬵�ǰ��Ԫ: %d', ie  ) ) ;
        k = StiffnessMatrix( ie ) ;
        AssembleStiffnessMatrix( ie, k ) ;
    end

    % step3. ����ֲ�ѹ�������ĵ�Ч�ڵ���
    for idfstep = 1:length(gDFStep)
        df_number = gDFStep( idfstep ) ;
        for i=1+sum(gDFStep(1:idfstep-1)):df_number+sum(gDFStep(1:idfstep-1))
            enf = EquivalentDistPressure( gDF(i,1), gDF(i,2), gDF(i,3), gDF(i,4) ) ;
            m = gDF(i, 1) ;
            n = gDF(i, 2) ;
            f(m*2-1:m*2,idfstep) = f(m*2-1:m*2,idfstep)+enf( 1:2 ) ;
            f(n*2-1:n*2,idfstep) = f(n*2-1:n*2,idfstep)+enf( 3:4 ) ;
        end
    end
    
    % step4. ����Լ���������޸ĸնȾ���ͽڵ������������û��л��з�
    [bc_number,dummy] = size( gBC1 ) ;
    for ibc=1:1:bc_number
        n = gBC1(ibc, 1 ) ;
        d = gBC1(ibc, 2 ) ;
        m = (n-1)*2 + d ;
        for idfstep = 1:length(gDFStep)
            f(:,idfstep) = f(:,idfstep) - gBC1(ibc,3) * gK(:,m) ;
            f(m,idfstep) = gBC1(ibc,3) ;
        end
        gK(:,m) = zeros( node_number*2, 1 ) ;
        gK(m,:) = zeros( 1, node_number*2 ) ;
        gK(m,m) = 1.0 ;
    end

    % step 5. ��ⷽ���飬�õ��ڵ�λ������
    gDelta = gK \ f ;
return

function k = StiffnessMatrix( ie )
%  ���㵥Ԫ�նȾ���
%  �������:
%     ie ----  ��Ԫ��
%  ����ֵ:
%     k  ----  ��Ԫ�նȾ���

    global gNode gElement gMaterial 
    E  = gMaterial( gElement(ie, 4), 1 ) ;
    mu = gMaterial( gElement(ie, 4), 2 ) ;
    ri = gNode( gElement( ie, 1 ), 1 ) ;
    zi = gNode( gElement( ie, 1 ), 2 ) ;
    rj = gNode( gElement( ie, 2 ), 1 ) ;
    zj = gNode( gElement( ie, 2 ), 2 ) ;
    rm = gNode( gElement( ie, 3 ), 1 ) ;
    zm = gNode( gElement( ie, 3 ), 2 ) ;
    ai = rj*zm - rm*zj ;
    aj = rm*zi - ri*zm ;
    am = ri*zj - rj*zi ;
    bi = zj - zm ;
    bj = zm - zi ;
    bm = zi - zj ;
    ci = -(rj-rm) ;
    cj = -(rm-ri) ;
    cm = -(ri-rj) ;
    area = (ai+aj+am)/2 ;
    r = (ri+rj+rm)/3; 
    z = (zi+zj+zm)/3 ;
    gi = ai/r + bi + ci*z/r ;
    gj = aj/r + bj + cj*z/r ;
    gm = am/r + bm + cm*z/r ;
    B = [bi  0 bj  0 bm  0
         gi  0 gj  0 gm  0
          0 ci  0 cj  0 cm
         ci bi cj bj cm bm] ;
    B = B/2/area ;
    A1 = mu/(1-mu) ;
    A2 = (1-2*mu)/2/(1-mu) ;
    A3 = E*(1-mu)/(1+mu)/(1-2*mu) ;
    D = A3*[ 1  A1  A1    0
            A1   1  A1    0
            A1  A1   1    0  
             0   0   0   A2] ;
    k = 2*pi*r*abs(area)*transpose(B)*D*B ;    
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

function enf = EquivalentDistPressure( m, n, pm, pn )
% �������Էֲ�ѹ���ĵ�Ч�ڵ���
% �������:
%      i,j ----------  ����
%      pi,pj --------  �����Ŷ�Ӧ��ѹ��ֵ
% ����ֵ:
%      enf -----------  ��Ч�ڵ�������
    global gNode
    enf = zeros( 4, 1 ) ;
    
    % ��ȡ�������
    rm = gNode( m, 1 ) ;
    rn = gNode( n, 1 ) ;
    zm = gNode( m, 2 ) ;
    zn = gNode( n, 2 ) ;

    % ����pi_, pj_
    pm_ = pi/6*( pm*(3*rm+rn) + pn*(rm+rn) ) ;
    pn_ = pi/6*( pm*(rm+rn) + pn*(rm+3*rn) ); 
    
    % ������������ϵĵ�Ч�ڵ���
    enf(1:2) = pm_ * [zm-zn;rn-rm] ;
    enf(3:4) = pn_ * [zm-zn;rn-rm] ;
return

function SaveResults( file_out )
%  ���������
%  ���������
%     file_out  --- �������ļ�
%  ����ֵ��
%     ��
    global gNode gElement gMaterial gBC1 gDelta
    save( file_out, 'gNode', 'gElement', 'gMaterial', 'gBC1', 'gDelta' ) ;
return

function DisplayResults
%  ��ʾ������
%  ���������
%     ��
%  ����ֵ��
%     ��
    global gNode gElement gMaterial gBC1 gDelta
    
    %  ������
    a = [ 1, 0.875, 0.750, 0.625, 0.5 ] ;
    width = 10 ;
    height = 15 ;
    q = 1e5 ;
    E = gMaterial( 1, 1 ) ;
    mu = gMaterial( 1, 2 ) ;
    w = -2*q*a*(1-mu^2)/E ;
    
    % ����Ԫ��
    node = 47 ;   % ԭ�㴦�Ľ���
    wfem = full(gDelta( node*2, : )) ;

    % ��ʾ���
    fprintf( '\n\n\n               ��ͬ�غ����ð뾶���غ�Բ�Ĵ�������λ��\n' ) ;
    fprintf( '=======================================================================\n' ) ;
    fprintf( '�غɰ뾶  ��������(m)      ����Ԫ        ������        ������\n' ) ;
    fprintf( '  (m)    (��ȡ����)       (mm)          (mm)             (%%) \n' ) ;
    fprintf( '-----------------------------------------------------------------------\n' ) ;
    for i=1:length(a)
        fprintf( '%6.3f     %2.0f��%2.0f       %8.3f       %8.3f        %8.3f\n', ...
                  a(i), height, width, wfem(i)*1000, w(i)*1000, (wfem(i)-w(i))/w(i)*100 ) ;
    end
    fprintf( '=======================================================================\n' ) ;
return

