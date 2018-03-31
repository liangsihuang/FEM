function exam7_1
% ������Ϊ�����µĵ�һ�����������ñ�����ε�Ԫ������ΰ������
%  ��������� 
%      ��

    % �����ĳߴ�
    a = 1 ;     % ��x������
    b = 1 ;     % ��y���򳤶�
    h = 0.01 ;  % ��z������
    
    % ���������ܶ�     
    m = 16 ;     % x����Ԫ��Ŀ 
    n = 16 ;     % y����Ԫ��Ŀ

    % �����������
    E = 2.1e11 ;   % ����ģ��
    mu = 0.3 ;     % ���ɱ�
    
    % ���ز���
    load = '������' ;   % Ҳ������ ������
    %load = '������' ;
    p = 1e6 ;             % �غɴ�С
    
    % �߽�����
    %bc = '��֧' ;  % Ҳ������ ��֧    
    bc = '��֧' ;
    
    % ����ļ�����
    file_out = 'exam7_1.mat' ;
    
    % �������ļ��Ƿ����
    file_prompt = 0 ;
    if exist( file_out ) ~= 0 & file_prompt == 1
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

    FemModel(a,b,h,m,n,E,mu,load,p,bc) ;       % ��������Ԫģ��
    SolveModel ;                     % �������Ԫģ��
    SaveResults( file_out ) ;        % ���������
    DisplayResults(a,b,h,m,n,E,mu,load,p,bc) ; % �Ѽ�������ʾ���� 
return

function FemModel(a,b,h,m,n,E,mu,load,p,bc)
%  ��������Ԫģ��
%  ���������
%    a  ---------   ��x������
%    b  ---------   ��y���򳤶�
%    h  ---------   ��z������
%    m  ---------   x����Ԫ��Ŀ 
%    n  ---------   y����Ԫ��Ŀ
%    E  ---------   ����ģ��
%    mu ---------   ���ɱ�
%    load -------   �غ�����
%    p  ---------   �غɴ�С
%    bc ---------   �߽�����
%  ����ֵ��
%      ��

    global gNode gElement gMaterial gBC1 gDF gNF gK

    % ����������
    dx = a / m / 2 ;   % ȡ���1/4���з���
    dy = b / n / 2 ;
    gNode = zeros( (m+1)*(n+1), 2 ) ;
    for i=1:n+1
        for j=1:m+1
            gNode( (i-1)*(m+1)+j, : ) = [dx*(j-1),dy*(i-1)] ;
        end
    end
    
    % ���嵥Ԫ
    gElement = zeros( n*n, 5 ) ;
    for i=1:n
        for j=1:m
            gElement( (i-1)*m+j, 1:4) = [ (i-1)*(m+1)+j, ...
                                          (i-1)*(m+1)+j+1, ...
                                          i*(m+1)+j+1,...
                                          i*(m+1)+j ] ;
        end
    end
    gElement( :, 5 ) = 1 ;
    
    % �������
    gMaterial = [ E, mu, h] ;  
    
    % ȷ���߽�����
    if bc == '��֧'
        gBC1 = zeros( (n+1)*2+m*2+1, 3 ) ;
        for i=1:m+1
            gBC1( i, : ) = [i, 2, 0] ;                    % x=0�ı߽�����x���ת�ǵ�����
        end
        for i=1:n+1
            gBC1( (m+1)+i, : ) = [(i-1)*(m+1)+1, 3, 0] ;  % y=0�ı߽�����y���ת�ǵ�����
        end
        for i=1:n+1
            gBC1( m+n+2+i, : ) = [i*(m+1), 1, 0 ] ;       % x=a/2�ı߽����Ӷ�Ϊ��
        end
        for i=1:m
            gBC1( (n+1)*2+m+1+i, : ) = [n*(m+1)+i, 1, 0] ;    % y=b/2�ı߽����Ӷ�Ϊ��
        end
    elseif bc == '��֧'
        gBC1 = zeros( m*4+(n+1)*3+n, 3 ) ;
        for i=1:m
            gBC1( i, : ) = [i, 2, 0] ;                    % x=0�ı߽�����x���ת�ǵ�����
        end
        for i=1:n
            gBC1( m+i, : ) = [(i-1)*(m+1)+1, 3, 0] ;      % y=0�ı߽�����y���ת�ǵ�����
        end
        for i=1:n+1
            gBC1( m+n+i, : ) = [i*(m+1), 1, 0 ] ;          % x=a/2�ı߽����Ӷ�Ϊ��
        end
        for i=1:n+1
            gBC1( m+n+n+1+i, : ) = [i*(m+1), 2, 0 ] ;      % x=a/2�ı߽�����x���ת��Ϊ��
        end
        for i=1:n+1
            gBC1( m+n+(n+1)*2+i, : ) = [i*(m+1), 3, 0 ] ;  % x=a/2�ı߽�����y���ת��Ϊ��
        end
        for i=1:m
            gBC1( m+n+(n+1)*3+i, : ) = [n*(m+1)+i, 1, 0] ; % y=b/2�ı߽����Ӷ�Ϊ��
        end
        for i=1:m
            gBC1( m*2+n+(n+1)*3+i, : ) = [n*(m+1)+i, 2, 0] ; % y=b/2�ı߽�����x���ת��Ϊ��
        end
        for i=1:m
            gBC1( m*3+n+(n+1)*3+i, : ) = [n*(m+1)+i, 3, 0] ; % y=b/2�ı߽�����y���ת��Ϊ��
        end
    end

    % ȷ���غ�
    if load == '������'
        gNF = [] ;
        gDF = gElement ;
        gDF(:,5) = p ;
    elseif load == '������'
        gDF = [] ;
        gNF = [1, 1, p] ;
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

    global gNode gElement gMaterial gBC1 gDF gNF gK gDelta gElementMoment gNodeMoment

    % step1. ��������նȾ���ͽڵ�������
    [node_number,dummy] = size( gNode ) ;
    gK = sparse( node_number * 3, node_number * 3 ) ;
    f = sparse( node_number * 3, 1) ;    

    % step2. ���㵥Ԫ�նȾ��󣬲����ɵ�����նȾ�����
    [element_number,dummy] = size( gElement ) ;
    for ie=1:1:element_number
        disp( sprintf(  '���㵥Ԫ�նȾ��󲢼��ɣ���ǰ��Ԫ: %d', ie  ) ) ;
        k = StiffnessMatrix( ie ) ;
        AssembleStiffnessMatrix( ie, k ) ;
    end
    
    % step3. ����ֲ�ѹ���ĵ�Ч�ڵ����������ɵ�����ڵ���������
    [df_number, dummy] = size( gDF ) ;
    for idf = 1:1:df_number
        edf = EquivalentDistPressure( gDF(idf,1:4), gDF(idf,5) ) ;
        node = gDF( idf, 1:4 ) ;
        f( node*3-2 ) = f( node*3-2 ) + edf( 1:3:10 ) ;
        f( node*3-1 ) = f( node*3-1 ) + edf( 2:3:11 ) ;
        f( node*3 )   = f( node*3 )   + edf( 3:3:12 ) ;
    end
  
    % step4. ����ֲ�ѹ���ĵ�Ч�ڵ����������ɵ�����ڵ���������
    [nf_number,dummy] = size( gNF ) ;
    for inf = 1:nf_number
        node = gNF( inf, 1 ) ;
        dim = gNF( inf, 2 ) ;
        f( (node-1)*3+dim ) = f( (node-1)*3+dim ) + gNF( inf, 3 ) ;
    end
    % step5. �����һ��Լ���������޸ĸնȾ���ͽڵ������������ó˴�����
    [bc1_number,dummy] = size( gBC1 ) ;
    for ibc=1:1:bc1_number
        n = gBC1(ibc, 1 ) ;
        d = gBC1(ibc, 2 ) ;
        m = (n-1)*3 + d ;
        f(m) = gBC1(ibc, 3)* gK(m,m) * 1e10 ;
        gK(m,m) = gK(m,m) * 1e10 ;
    end
    
    % step5. ��ⷽ���飬�õ��ڵ�λ������
    gDelta = gK \ f ;
    
    % step6. ���㵥Ԫ���
    gElementMoment = zeros( element_number, 4, 3 ) ;
    for ie=1:element_number 
        disp( sprintf(  '���㵥Ԫ��أ���ǰ��Ԫ: %d', ie  ) ) ;
        em = ElementMoment( ie ) ;
        gElementMoment( ie, :, : ) = em ;
    end
    
    % step7. ����ڵ����(�����ƽڵ�ƽ��)
    gNodeMoment = zeros( node_number, 3 ) ;       
    nodenum = zeros( node_number,1 ) ;
    for i=1:node_number                         
        disp( sprintf(  '����ڵ���أ���ǰ���: %d', i  ) ) ;
        for ie=1:element_number
            index = find(  gElement( ie, 1:4 ) == i ) ;
            if ~isempty(index)
                nodenum(i) = nodenum(i) + 1 ;
                gNodeMoment(i,1) = gNodeMoment(i,1) + gElementMoment( ie, index, 1 ) ;
                gNodeMoment(i,2) = gNodeMoment(i,2) + gElementMoment( ie, index, 2 ) ;
                gNodeMoment(i,3) = gNodeMoment(i,3) + gElementMoment( ie, index, 3 ) ;
            end
        end
    end
    gNodeMoment(:,1) = gNodeMoment(:,1) ./ nodenum ;
    gNodeMoment(:,2) = gNodeMoment(:,2) ./ nodenum ;
    gNodeMoment(:,3) = gNodeMoment(:,3) ./ nodenum ;
return

function k = StiffnessMatrix( ie )   
%  ����ƽ��Ӧ��Ȳ�����Ԫ�ĸնȾ���
%  ���������
%     ie -- ��Ԫ��
%  ����ֵ��
%     k  -- ��Ԫ�նȾ���

    global gNode gElement gMaterial
    k = zeros( 12, 12 ) ;
    E = gMaterial( gElement( ie, 5 ), 1 ) ;
    mu = gMaterial( gElement( ie, 5 ), 2 ) ;
    h = gMaterial( gElement( ie, 5 ), 3 ) ;
    D = E * h^3 / 12 / ( 1 - mu^2 ) ;
    x = gNode( gElement( ie, : ), 1 ) ;
    y = gNode( gElement( ie, : ), 2 ) ;
    a = abs( x(3) - x(1) ) / 2 ;
    b = abs( y(3) - y(1) ) / 2 ;
    ab2 = a^2/b^2 ;
    ba2 = b^2/a^2 ;
    H = D/60/a/b ;
    xi = [ -1, 1, 1, -1 ] ;
    eta = [ -1, -1, 1, 1 ] ;
    for i=1:4
        for j=1:4
            xi0 = xi(i)*xi(j) ;
            eta0 = eta(i)*eta(j) ;
            a11 = 3*H*(15*(ba2*xi0+ab2*eta0) ...
                   +(14-4*mu+5*ba2+5*ab2)*xi0*eta0) ;
            a12 = -3*H*b*((2+3*mu+5*ab2)*xi0*eta(i) ...
                   + 15*ab2*eta(i)+5*mu*xi0*eta(j)) ;
            a13 = 3*H*a*((2+3*mu+5*ba2)*xi(i)*eta0 ...
                   + 15*ba2*xi(i)+5*mu*xi(j)*eta0) ;
            a21 = -3*H*b*((2+3*mu+5*ab2)*xi0*eta(j) ...
                   + 15*ab2*eta(j)+5*mu*xi0*eta(i)) ;
            a22 = H*b^2*(2*(1-mu)*xi0*(3+5*eta0)+...
                   + 5*ab2*(3+xi0)*(3+eta0)) ;
            a23 = -15*H*mu*a*b*(xi(i)+xi(j))*(eta(i)+eta(j)) ;
            a31 = 3*H*a*((2+3*mu+5*ba2)*xi(j)*eta0 ...
                   + 15*ba2*xi(j)+5*mu*xi(i)*eta0) ;
            a32 = -15*H*mu*a*b*(xi(i)+xi(j))*(eta(i)+eta(j)) ;
            a33 = H*a^2*(2*(1-mu)*eta0*(3+5*xi0)+...
                   + 5*ba2*(3+xi0)*(3+eta0)) ;
            k( i*3-2:i*3, j*3-2:j*3 ) = [ a11 a12 a13;
                                          a21 a22 a23;
                                          a31 a32 a33 ] ;
        end
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

function edf = EquivalentDistPressure( node, press )
%   ����ֲ�ѹ���ĵ�Ч�ڵ���
%   �������:
%      node ----------  ����
%      press ---------  �����Ŷ�Ӧ��ѹ��ֵ
%   ����ֵ:
%      edf -----------  ��Ч�ڵ�������
    global gNode
    
    edf = zeros( 12, 1 ) ;
    x = gNode( node, 1 ) ;
    y = gNode( node, 2 ) ;
    a = abs( x(3)-x(1) ) / 2 ;
    b = abs( y(3)-y(1) ) / 2 ;
    xi = [-1; 1; 1; -1] ;
    eta = [-1;-1; 1; 1] ;
    for i=1:4
        edf(i*3-2:i*3) = [press*a*b; -press*a*b^2/3*eta(i); press*a^2*b/3*xi(i)] ;
    end
return

function em = ElementMoment( ie )
%  ���㵥Ԫ�Ľ�����
%  �������:
%      ie ----- ��Ԫ��
%  ����ֵ:
%      em ----- ��Ԫ������

    global gElement gDelta gMaterial

    E = gMaterial( gElement( ie, 5 ),  1 ) ;
    mu = gMaterial( gElement( ie, 5 ), 2 ) ;
    h = gMaterial( gElement( ie, 5 ),  3 ) ;
    em = zeros( 4, 3 ) ;
    node = gElement( ie, 1:4 ) ;
    delta = zeros( 12, 1 ) ;
    delta(1:3:10) = gDelta( (node-1)*3+1 ) ;
    delta(2:3:11) = gDelta( (node-1)*3+2 ) ;
    delta(3:3:12) = gDelta( (node-1)*3+3 ) ;
    xi  = [ -1, 1, 1, -1 ] ;
    eta = [ -1, -1, 1, 1 ] ;
    D = E/(1-mu^2)*[1  mu  0
                    mu  1  0
                    0   0  (1-mu)/2 ] ;
    for i=1:4
        B = MatrixB( ie, xi(i), eta(i) ) ;
        em( i, : ) = transpose(h^3/12*D*B*delta) ;
    end
return

function B = MatrixB( ie, x, e )
%  ���㵥Ԫ��Ӧ�����B
%  �������:
%      ie ----- ��Ԫ��
%  ����ֵ:
%      B ------ ��ԪӦ�����
    
    global gNode gElement
    B = zeros( 3, 12 ) ;
    xy = gNode( gElement(ie, 1:4), : ) ;
    a = abs(xy(3,1)-xy(1,1))/2 ;
    b = abs(xy(3,2)-xy(1,2))/2 ;
    xi  = [ -1, 1, 1, -1 ] ;
    eta = [ -1, -1, 1, 1 ] ;
    for i=1:4
        x0 = xi(i)*x ;
        e0 = eta(i)*e ;
        B(:,(i-1)*3+1:(i-1)*3+3) = [ 3*b/a*x0*(1+e0)                       0                      b*xi(i)*(1+3*x0)*(1+e0)
                                     3*a/b*e0*(1+x0)                -a*eta(i)*(1+x0)*(1+3*e0)           0
                                     xi(i)*eta(i)*(3*x^2+3*e^2-4)   -b*xi(i)*(3*e^2+2*e0-1)       a*eta(i)*(3*x^2+2*x0-1) ] ;
    end
    B = B / a / b / 4 ;
return

function SaveResults( file_out )
%  ���������
%  ���������
%     ��
%  ����ֵ��
%     ��
    global gNode gElement gMaterial gBC1 gDF gDelta gNF gK
    save( file_out, 'gNode', 'gElement', 'gMaterial', 'gBC1',  ...
          'gDF', 'gDelta', 'gNF', 'gK' ) ;
return

function DisplayResults(a,b,h,m,n,E,mu,load,p,bc)
%  ��ʾ��������������������Ƚ�
%  ���������
%    a  ---------   ��x������
%    b  ---------   ��y���򳤶�
%    h  ---------   ��z������
%    m  ---------   x����Ԫ��Ŀ 
%    n  ---------   y����Ԫ��Ŀ
%    E  ---------   ����ģ��
%    mu ---------   ���ɱ�
%    load -------   �غ�����
%    p  ---------   �غɴ�С
%    bc ---------   �߽�����
%  ����ֵ��
%     ��

    global gDelta gNodeMoment
    
    wmax = gDelta( 1, 1 ) ;          % �����Ӷ�
    D = E*h^3/12/(1-mu^2) ;          % ����ն�

    mx  = gNodeMoment( 1, 1 ) ;   % ����������(Mx)
    my  = gNodeMoment( 1, 2 ) ;   % ����������(My)
    mxy = gNodeMoment( 1, 3 ) ;   % ����������(Mxy)
    
    fprintf( '��ĳߴ�a��b��h = %g��%g��%g\n', a, b, h ) ;
    fprintf( '����ģ�� = %g\n', E ) ;
    fprintf( '���ɱ�   = %g\n', mu ) ;
    fprintf( '�غ����� = %s\n', load ) ;
    fprintf( '�غɴ�С = %g\n', p ) ;
    fprintf( '�߽����� = %s\n', bc ) ;
    fprintf( '�����ܶ� m��n = %d��%d\n', m, n ) ;
    if load == '������'
        alpha = wmax * D / p / a^4 ;
        beta = mx/p/a^2 ;
        beta1 = my/p/a^2 ;
    elseif load == '������'
        alpha = wmax * D / p / 4 / a^2 ;
        beta = mx/p/4 ;
        beta1 = my/p/3 ;
    end
    
    fprintf( '�Ӷ�ϵ�� alpha = %.7f\n', alpha ) ;
    fprintf( '���ϵ�� beta  = %.7f\n', beta ) ;
    fprintf( '���ϵ�� beta1 = %.7f\n', beta1 ) ;
%    fid = fopen( 'exam7_1.out', 'a' ) ;
%        fprintf( fid, '%.1f  %15.7f  %15.7f\n', b/a, alpha, beta ) ;
%    fclose( fid ) ;
return