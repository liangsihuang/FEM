function exam4_1( m, n )
% ������Ϊ�����µĵ�һ�����������þ��ε�Ԫ���㴿�����ı���
%      exam4_1(m,n) 
%  ��������� 
%      m ------  x����Ԫ��Ŀ
%      n ------  y����Ԫ��Ŀ

% ����ȫ�ֱ���
%      gNode ------ �ڵ�����
%      gElement --- ��Ԫ����
%      gMaterial -- ��������
%      gBC1 ------- ��һ��Լ������
%      gDF -------- �ֲ���
%      gK --------- ����նȾ���
%      gDelta ----- ����ڵ�����
    global gNode gElement gMaterial gBC1 gDF gK gDelta

    if nargin < 1 
        m = 4 ;
        n = 4 ;
    elseif nargin < 2
        n = 4 ;
    end
    
    FemModel(m, n) ;       % ��������Ԫģ��
    SolveModel ;           % �������Ԫģ��
    DisplayResults ;       % ��ʾ������
return ;

function FemModel(m, n)
%  ��������Ԫģ��
%  ���������
%      m ---  x����Ԫ��Ŀ
%      n ---  y����Ԫ��Ŀ
%  ����ֵ��
%      ��
%  ˵����
%      �ú�������ƽ���ϵ������Ԫģ�����ݣ�
%        gNode ------- �ڵ㶨��
%        gElement ---- ��Ԫ����
%        gMaterial --- ���϶��壬��������ģ�������Ľ���������Ŀ�����Ծ�
%        gBC --------- Լ������
%        gDF --------- �ֲ���

    global gNode gElement gMaterial gBC1 gDF

    length = 0.045 ;    % ���㲿�ֵĳ���(x����)
    height = 0.03 ;     % ���㲿�ֵĸ߶�(y����)
    dx = length/m ;     % ���ε�Ԫ�Ŀ��
    dy = height/n ;     % ���ε�Ԫ�ĸ߶�
    p = 10e6 ;          % ����Ӧ��10MPa
    
    % �ڵ�����
    gNode = zeros( (m+1)*(n+1), 2 ) ;
    for i=1:n+1
        for j=1:m+1
            k = (i-1)*(m+1)+j ;        % �ڵ��
            xk = (j-1)*dx ;            % �ڵ��x����
            yk = (i-1)*dy ;            % �ڵ��y����
            gNode(k,:) = [xk, yk] ;
        end
    end
     
    % ��Ԫ����
    gElement = zeros( m*n, 4 ) ;
    for i=1:n
        for j=1:m
            k = (i-1)*m+j ;             % ��Ԫ��
            n1 = (i-1)*(m+1)+j ;        % ��һ���ڵ��
            n2 = (i-1)*(m+1)+j+1 ;      % �ڶ����ڵ��
            n3 = i*(m+1)+j+1 ;          % �������ڵ��
            n4 = i*(m+1)+j ;            % ���ĸ��ڵ��
            gElement(k,:) = [n1, n2, n3, n4] ;
        end
    end

    % �������� 
    %           ����ģ��    ���ɱ�   ���
    gMaterial = [2.0e11,    0.3,     0.01] ;   %  ���� 1

    % ��һ��Լ������
    gBC1 = zeros( m+1+n+1, 3 ) ;
    for j=1:m+1
        gBC1(j,:) = [j, 1, 0.0] ;
    end
    for i=2:n+1
        gBC1(m+1+i,:) = [(i-1)*(m+1)+1, 1, 0.0] ; 
    end
    gBC1(m+1+1,:) = [1, 2, 0.0] ;

    % �ֲ��غɣ����Էֲ���
    gDF = zeros( n, 5 ) ;
    for i=1:n
        k = i*m ;
        gDF(i,:) = [ k, 2, (i-1)*p/n, i*p/n, 1] ;
    end
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

    global gNode gElement gMaterial gBC1 gDF gK gDelta

    % step1. ��������նȾ���ͽڵ�������
    [node_number,dummy] = size( gNode ) ;
    gK = sparse( node_number * 2, node_number * 2 ) ;
    f = sparse( node_number * 2, 1 ) ;

    % step2. ���㵥Ԫ�նȾ��󣬲����ɵ�����նȾ�����
    [element_number,dummy] = size( gElement ) ;
    for ie=1:1:element_number
        k = StiffnessMatrix( ie ) ;
        AssembleStiffnessMatrix( ie, k ) ;
    end

    % step3. ����ֲ����ĵ�Ч�ڵ����������ɵ�����ڵ���������
    [df_number, dummy] = size( gDF ) ;
    for idf = 1:1:df_number
        enf = EquivalentNodeForce( gDF(idf,1), gDF(idf,2), gDF(idf,3), gDF(idf,4), gDF(idf,5) ) ;
        ielem = gDF(idf,1) ;
        iedge = gDF(idf,2) ;
        i = gElement( ielem, iedge ) ;
        if iedge < 4 
            j = gElement( ielem, iedge+1 ) ;
        else
            j = gElement( ielem, 1 );
        end
        
        f( (i-1)*2+1 : (i-1)*2+2 ) = f( (i-1)*2+1 : (i-1)*2+2 ) + enf( 1:2 ) ;
        f( (j-1)*2+1 : (j-1)*2+2 ) = f( (j-1)*2+1 : (j-1)*2+2 ) + enf( 3:4 ) ;
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
return


function k = StiffnessMatrix( ie )
%  ���㵥Ԫ�նȾ���
%  �������:
%     ie ----  ��Ԫ��
%  ����ֵ:
%     k  ----  ��Ԫ�նȾ���

global gNode gElement gMaterial 
    k = zeros( 8, 8 ) ;
    E  = gMaterial( 1 ) ;
    mu = gMaterial( 2 ) ;
    h  = gMaterial( 3 ) ;
    x1 = gNode( gElement( ie, 1 ), 1 ) ;
    y1 = gNode( gElement( ie, 1 ), 2 ) ;
    x3 = gNode( gElement( ie, 3 ), 1 ) ;
    y3 = gNode( gElement( ie, 3 ), 2 ) ;
    a = (x3-x1)/2 ;
    b = (y3-y1)/2 ;
    xi  = [ -1, 1, 1, -1 ] ;
    eta = [ -1, -1, 1, 1 ] ;
    for i=1:1:4
        for j=1:1:4
            k( (i-1)*2 + 1, (j-1)*2 + 1 ) = b/a*xi(i)*xi(j)*(1+1/3*eta(i)*eta(j)) + ...
                (1-mu)/2*a/b*eta(i)*eta(j)*(1+1/3*xi(i)*xi(j)); 
            k( (i-1)*2 + 1, (j-1)*2 + 2 ) = mu*eta(j)*xi(i) + (1-mu)/2*xi(j)*eta(i) ;
            k( (i-1)*2 + 2, (j-1)*2 + 1 ) = mu*xi(j)*eta(i) + (1-mu)/2*eta(j)*xi(i) ;
            k( (i-1)*2 + 2, (j-1)*2 + 2 ) = a/b*eta(i)*eta(j)*(1+1/3*xi(i)*xi(j)) + ...
                (1-mu)/2*b/a*xi(i)*xi(j)*(1+1/3*eta(i)*eta(j)) ;
        end
    end
    eh = E*h/4/(1-mu^2) ;
    k = eh*k ;
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

function enf = EquivalentNodeForce( ie, iedge, p1, p2, idof )
%   �������Էֲ����صĵ�Ч�ڵ���
%   �������:
%      ie  -----  ��Ԫ��
%      ieged ---  ���õıߺ�
%      p1  -----  ��һ���ڵ��ϵķֲ�������ֵ
%      p2  -----  �ڶ����ڵ��ϵķֲ�������ֵ
%      idof  ---  �ֲ����ķ���
%                  1 ---  x����
%                  2 ---  y����
%   ����ֵ:
%      enf -----  ��Ч�ڵ�������
    global gElement gNode gMaterial
    h = gMaterial( 3 ) ;
    x1 = gNode( gElement( ie, 1 ), 1 ) ;
    y1 = gNode( gElement( ie, 1 ), 2 ) ;
    x3 = gNode( gElement( ie, 3 ), 1 ) ;
    y3 = gNode( gElement( ie, 3 ), 2 ) ;
    a = ( x3 - x1 ) / 2 ;
    b = ( y3 - y1 ) / 2 ;
    if iedge == 1 | iedge == 3
        length = a ;        
    else
        length = b ;
    end
    f1 = h*length*(2*p1+p2)/3 ;
    f2 = h*length*(p1+2*p2)/3 ;
    if idof == 1
        enf = [f1; 0; f2; 0] ;
    else
        enf = [0; f1; 0; f2] ;
    end
return


function DisplayResults
%  ��ʾ������
%  ���������
%     ��
%  ����ֵ��
%     ��

    global gNode gDelta
    
    fprintf( '�ڵ�λ��\n' ) ; 
    fprintf( '  �ڵ��         x����λ��               y����λ��\n' ) ; 
    [node_number,dummy] = size( gNode ) ;
    for i=1:node_number
        fprintf(  '%6d       %16.8e        %16.8e\n',...
                  i, gDelta((i-1)*2+1), gDelta((i-1)*2+2) ) ; 
    end
return
