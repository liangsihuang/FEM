function exam3_2
% ������Ϊ�����µĵڶ�������������ƽ������Ԫ����б�ȸռ��ŵı��κ�����
%   ��������� ��
%   �������� �ڵ�λ�ƺ͵�Ԫ�ڵ��� 

    % ��������ļ��Ƿ����
    file_in = 'exam3_2.dat' ;
    if exist( file_in ) == 0
        disp( sprintf( '�����ļ� %s ������', file_in ) )
        disp( sprintf( '������ֹ' ) )
        return ;
    end

    % ����ģ�ͣ���Ⲣ��ʾ���
    PlaneFrameModel(file_in) ;   % ��������Ԫģ��
    SolveModel ;                 % �������Ԫģ��
    DisplayResults ;             % ��ʾ������
return ;

function PlaneFrameModel( file_in )
%  ����ƽ���ϵ������Ԫģ��
%  ���������
%      file_in ------- ����Ԫģ�������ļ�
%  ����ֵ��
%      ��
%  ˵����
%      �ú�������ƽ���ϵ������Ԫģ�����ݣ�
%        gNode ------- �ڵ㶨��
%        gElement ---- ��Ԫ����
%        gMaterial --- ���϶��壬��������ģ�������Ľ���������Ŀ�����Ծ�
%        gBC1 -------- ��һ��߽�����
%        gBC2 -------- �ڶ���Լ������, ����½ӵĽ��

    global gNode gElement gMaterial gBC1 gBC2

    % �������ļ�
    fid_in = fopen( file_in, 'r' ) ;
    
    % ��ȡ�ڵ�����ͽڵ��¶�
    node_number = fscanf( fid_in, '%d', 1 ) ;
    gNode = zeros( node_number, 3 ) ;
    for i=1:node_number
        dummy = fscanf( fid_in, '%d', 1 ) ;
        gNode(i,:) = fscanf( fid_in, '%f', [1 3] ) ;
    end
     
    % ��ȡ��Ԫ����
    element_number = fscanf( fid_in, '%d', 1 ) ;
    gElement = zeros( element_number, 3 ) ;
    for i=1:element_number
        dummy = fscanf( fid_in, '%d', 1 ) ;
        gElement(i,:) = fscanf( fid_in, '%d', [1 3] ) ;
    end
        
    % ��ȡ�������� 
    material_number = fscanf( fid_in, '%d', 1 ) ;
    gMaterial = zeros( material_number, 5 ) ;
    for i=1:material_number
        dummy = fscanf( fid_in, '%d', 1 ) ;
        gMaterial(i,:) = fscanf( fid_in, '%f', [1 5] ) ;
    end
    
    % ��ȡ��һ��Լ������
    bc1_number = fscanf( fid_in, '%d', 1 ) ;
    gBC1 = zeros( bc1_number, 3 ) ;
    for i=1:bc1_number
        gBC1(i,1) = fscanf( fid_in, '%d', 1 ) ;
        gBC1(i,2) = fscanf( fid_in, '%d', 1 ) ;
        gBC1(i,3) = fscanf( fid_in, '%f', 1 ) ;
    end

    % ��ȡ�ڶ���Լ������
    bc2_number = fscanf( fid_in, '%d', 1 ) ;
    gBC2 = zeros( bc2_number, 3 ) ;
    for i=1:bc2_number
        gBC2(i,:) = fscanf( fid_in, '%d', [1 3] ) ;
    end
    
    % �ر������ļ�
    fclose( fid_in ) ;
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

    global gNode gElement gMaterial gBC1 gBC2 gK gDelta

    % step1. ��������նȾ���ͽڵ�������
    [node_number,dummy] = size( gNode ) ;
    gK = sparse( node_number * 3, node_number * 3 ) ;
    f = sparse( node_number * 3, 1 ) ;

    % step2. ���㵥Ԫ�նȾ��󣬲����ɵ�����նȾ�����
    [element_number,dummy] = size( gElement ) ;
    for ie=1:1:element_number
        k = StiffnessMatrix( ie, 1 ) ;
        AssembleStiffnessMatrix( ie, k ) ;
    end

    % step3. ���������¶ȱ仯�����ĵ�Ч�ڵ���
    for ie=1:1:element_number
        egf = EquivalentGravityForce( ie ) ;
        etf = EquivalentThermalForce( ie ) ;
        i = gElement( ie, 1 ) ;
        j = gElement( ie, 2 ) ;
        f( (i-1)*3+1 : (i-1)*3+3 ) = f( (i-1)*3+1 : (i-1)*3+3 ) + egf( 1:3 ) + etf( 1:3 );
        f( (j-1)*3+1 : (j-1)*3+3 ) = f( (j-1)*3+1 : (j-1)*3+3 ) + egf( 4:6 ) + etf( 4:6 );
    end
    
    % step4. ����ڶ���Լ���������޸ĸնȾ���ͽڵ��������������û��л��з���
    [bc2_number,dummy] = size( gBC2 ) ;
    for ibc=1:1:bc2_number
        n1 = gBC2(ibc, 1 ) ;
        n2 = gBC2(ibc, 2 ) ;
        d = gBC2(ibc, 3 ) ;
        m1 = (n1-1)*3 + d ;
        m2 = (n2-1)*3 + d ;
        
        f(m2) = f(m2)+f(m1) ;
        f(m1) = 0 ;
        
        gK(m2,:) = gK(m2,:)+gK(m1,:) ;
        gK(:,m2) = gK(:,m2)+gK(:,m1) ;
        
        gK(m1,:) = zeros(1, node_number*3) ;
        gK(:,m1) = zeros(node_number*3, 1) ;
        gK(m1,m1) = 1.0 ;
        gK(m1,m2) = -1 ;
        gK(m2,m1) = -1 ;
        gK(m2,m2) = gK(m2,m2)+1 ;
    end
    
    % step5. �����һ��Լ���������޸ĸնȾ���ͽڵ��������������û��л��з���
    [bc1_number,dummy] = size( gBC1 ) ;
    for ibc=1:1:bc1_number
        n = gBC1(ibc, 1 ) ;
        d = gBC1(ibc, 2 ) ;
        m = (n-1)*3 + d ;
        f = f - gBC1(ibc,3) * gK(:,m) ;
        f(m) = gBC1(ibc,3) ;
        gK(:,m) = zeros( node_number*3, 1 ) ;
        gK(m,:) = zeros( 1, node_number*3 ) ;
        gK(m,m) = 1.0 ;
    end

    % step 7. ��ⷽ���飬�õ��ڵ�λ������
    gDelta = gK \ f ;
return

function DisplayResults
%  ��ʾ������
%  ���������
%     ��
%  ����ֵ��
%     ��

    global gNode gElement gMaterial gBC1 gNF gDF gK gDelta
    
    fprintf( '�ڵ�λ��\n' ) ; 
    fprintf( '  �ڵ��         x����λ��               y����λ��               ת��\n' ) ; 
    [node_number,dummy] = size( gNode ) ;
    for i=1:node_number
        fprintf(  '%6d       %16.8e        %16.8e       %16.8e\n',...
                  i, gDelta((i-1)*3+1), gDelta((i-1)*3+2), gDelta((i-1)*3+3) ) ; 
    end
    fprintf( '\n\n�ڵ���\n' ) ; 
    fprintf( '                                        ����               ����               ���\n' ) ; 
    [element_number, dummy] = size( gElement ) ;
    for ie = 1:element_number
        enf = ElementNodeForce( ie ) ;
        fprintf( '��Ԫ��%6d    �ڵ��%6d     %16.8e  %16.8e  %16.8e\n', ...
                  ie, gElement(ie,1), enf(1), enf(2), enf(3) ) ;
        fprintf( '                �ڵ��%6d     %16.8e  %16.8e  %16.8e\n', ...
                  gElement(ie,2), enf(4), enf(5), enf(6) ) ;
    end
return

function k = StiffnessMatrix( ie, icoord )
%  ���㵥Ԫ�նȾ���
%  �������:
%     ie -------  ��Ԫ��
%     icoord  --  ����ϵ��������������������֮һ
%                    1  ----  ��������ϵ
%                    2  ----  �ֲ�����ϵ
%  ����ֵ:
%     k  ----  ����icoord��ֵ����Ӧ����ϵ�µĸնȾ���
    global gNode gElement gMaterial
    k = zeros( 6, 6 ) ;
    E = gMaterial( gElement(ie, 3), 1 ) ;
    I = gMaterial( gElement(ie, 3), 2 ) ;
    A = gMaterial( gElement(ie, 3), 3 ) ;
    xi = gNode( gElement( ie, 1 ), 1 ) ;
    yi = gNode( gElement( ie, 1 ), 2 ) ;
    xj = gNode( gElement( ie, 2 ), 1 ) ;
    yj = gNode( gElement( ie, 2 ), 2 ) ;
    L = ( (xj-xi)^2 + (yj-yi)^2 )^(1/2) ;
    k = [  E*A/L           0          0 -E*A/L           0          0
               0  12*E*I/L^3  6*E*I/L^2      0 -12*E*I/L^3  6*E*I/L^2
               0   6*E*I/L^2    4*E*I/L      0  -6*E*I/L^2    2*E*I/L
          -E*A/L           0          0  E*A/L           0          0
               0 -12*E*I/L^3 -6*E*I/L^2      0  12*E*I/L^3 -6*E*I/L^2
               0   6*E*I/L^2    2*E*I/L      0  -6*E*I/L^2    4*E*I/L] ;
    if icoord == 1
        T = TransformMatrix( ie ) ;
        k = T*k*transpose(T) ;
    end
return

function AssembleStiffnessMatrix( ie, k )
%  �ѵ�Ԫ�նȾ��󼯳ɵ�����նȾ���
%  �������:
%      ie  --- ��Ԫ��
%      k   --- ��Ԫ�նȾ���
%  ����ֵ:
%      ��
    global gElement gK
    for i=1:1:2
        for j=1:1:2
            for p=1:1:3
                for q =1:1:3
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

function T = TransformMatrix( ie )
%  ���㵥Ԫ������ת������( �ֲ����� -> �������� )
%  �������
%      ie  ----- �ڵ��
%  ����ֵ
%      T ------- �Ӿֲ����굽�������������ת������
    global gElement gNode
    xi = gNode( gElement( ie, 1 ), 1 ) ;
    yi = gNode( gElement( ie, 1 ), 2 ) ;
    xj = gNode( gElement( ie, 2 ), 1 ) ;
    yj = gNode( gElement( ie, 2 ), 2 ) ;
    L = sqrt( (xj-xi)^2 + (yj-yi)^2 ) ;
    c = (xj-xi)/L ;
    s = (yj-yi)/L ;
    T=[ c  -s   0   0   0   0
        s   c   0   0   0   0
        0   0   1   0   0   0
        0   0   0   c  -s   0
        0   0   0   s   c   0
        0   0   0   0   0   1] ;
return

function egf = EquivalentGravityForce( ie )
%  ���㵥Ԫ���صĵ�Ч�ڵ���
%  �������
%      ie  ----- �ڵ��
%  ����ֵ
%      egf ----- ��������ϵ�µĵ�Ч�ڵ��� 
    global gElement gNode gMaterial
    
    xi = gNode( gElement( ie, 1 ), 1 ) ;
    yi = gNode( gElement( ie, 1 ), 2 ) ;
    xj = gNode( gElement( ie, 2 ), 1 ) ;
    yj = gNode( gElement( ie, 2 ), 2 ) ;
    L = sqrt( (xj-xi)^2 + (yj-yi)^2 ) ;
    c = (xj-xi)/L ;
    
    A  = gMaterial( gElement( ie, 3 ), 3 ) ;
    ro = gMaterial( gElement( ie, 3 ), 4 ) ;
    g = 9.8 ;
    
    egf = [ 0; -ro*g*A*L/2; -ro*g*A*L^2*c/12; 0; -ro*g*A*L/2; ro*g*A*L^2*c/12 ] ;
return

function etf = EquivalentThermalForce( ie )
%  ���㵥Ԫ���¶Ⱥ��صĵ�Ч�ڵ���
%  �������
%      ie  ----- �ڵ��
%  ����ֵ
%      etf ----- ��������ϵ�µĵ�Ч�ڵ��� 
    global gElement gNode gMaterial
    
    dT1 = gNode( gElement( ie, 1 ), 3 ) ;
    dT2 = gNode( gElement( ie, 2 ), 3 ) ;
    E = gMaterial( gElement( ie, 3 ), 1 ) ;
    A = gMaterial( gElement( ie, 3 ), 3 ) ;
    alpha = gMaterial( gElement( ie, 3 ), 5 ) ;
    Nx = E*A*alpha ;
    
    etf = [ -Nx; 0; 0; Nx; 0; 0 ] ;
    T = TransformMatrix( ie ) ;
    etf = T * etf ;
    etf = (dT1+dT2)/2 * etf ;
return

function enf = ElementNodeForce( ie )
%  ���㵥Ԫ�Ľڵ���
%  �������
%      ie  ----- �ڵ��
%  ����ֵ
%      enf ----- ��Ԫ�ֲ�����ϵ�µĽڵ��� 
    global gElement gNode gDelta
    i = gElement( ie, 1 ) ;
    j = gElement( ie, 2 ) ;
    de = zeros( 6, 1 ) ;
    de( 1:3 ) = gDelta( (i-1)*3+1:(i-1)*3+3 ) ;
    de( 4:6 ) = gDelta( (j-1)*3+1:(j-1)*3+3 ) ;
    k = StiffnessMatrix( ie, 1 ) ;
    enf = k * de ;
    etf = EquivalentThermalForce( ie ) ;
    egf = EquivalentGravityForce( ie ) ;
    enf = enf - etf -egf ;
    T = TransformMatrix( ie ) ;
    enf = transpose( T ) * enf ;
return