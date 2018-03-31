function exam6_1
% ������Ϊ�����µĵ�һ������������4��8����ı��εȲ�Ԫ��������ѹ����ת���ԲͲ
%  ��������� 
%       ��

    % ����ļ�����
    file_out = 'exam6_1.mat' ;
    
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

    FemModel ;                  % ��������Ԫģ��
    SolveModel ;                % �������Ԫģ��
    SaveResults( file_out ) ;   % ���������
    DisplayResults ;            % �Ѽ�������ʾ���� 
return

function FemModel
%  ��������Ԫģ��
%  ���������
%      ��
%  ����ֵ��
%      ��

%  ȫ�ֱ���������
%      gNode ------ �ڵ�����
%      gElement --- ��Ԫ����
%      gMaterial -- ��������
%      gBC1 ------- ��һ��Լ������
%      gBC3 ------- бԼ��
%      gDF -------- �ֲ���
    global gNode gElement gMaterial gBC1 gBC3 gDF
    
    r1 = 0.05 ;              %  �ڰ뾶
    r2 = 0.1 ;               %  ��뾶
    theta = 10.0*pi/180.0 ;  %  Բ�Ľǣ����ȣ�
    
    % �ڵ�����
    node_number = 28 ;
    gNode = zeros( node_number, 2 ) ;
    j = 0 ;
    for i = 1:5:26     % ���1,6,11,16,21,26������
        r = r1 + (r2-r1)/5 * j ;
        gNode(i,1) = r*cos(theta) ;
        gNode(i,2) = r*sin(theta) ;
        j = j + 1 ;
    end
    j = 0 ;
    for i = 4:5:24     % ���4,9,14,19,24������
        r = r1+(r2-r1)/10+(r2-r1)/5*j ;
        gNode(i,1) = r*cos(theta) ;
        gNode(i,2) = r*sin(theta) ;
        j = j+1 ;
    end
    j = 0 ;
    for i=2:5:27       % ���2,7,12,17,22,27������
        r = r1+(r2-r1)/5*j ;
        gNode(i,1) = r*cos(theta/2) ;
        gNode(i,2) = r*sin(theta/2) ;
        j = j+1 ;
    end
    j = 0 ;
    for i=3:5:28      % ���3,8,13,18,23,28������
        r = r1+(r2-r1)/5*j ;
        gNode(i,1) = r ;
        gNode(i,2) = 0 ;
        j = j+1 ;
    end
    j = 0 ;
    for i=5:5:25      % ���5,10,15,20,25������
        r = r1+(r2-r1)/10+(r2-r1)/5*j ;
        gNode(i,1) = r ;
        gNode(i,2) = 0 ;
        j = j+1 ;
    end
    
    % ��Ԫ����
    element_number = 5 ;
    gElement = zeros( element_number, 9 ) ;
    gElement(1,:)=[3  8  6  1   5  7  4  2  1] ;
    for i=2:5
        gElement(i,1:8) = gElement(i-1,1:8)+5 ;
        gElement(i,9) = 1 ;
    end
    
    % ��������
    gMaterial = [2.0e11 0.3 7800] ;
    
    % ��һ��Լ������
    bc1_number = 11;
    gBC1 = zeros( bc1_number, 3 ) ;
    j = 0 ;
    for i=3:5:28      %  ���3,8,13,18,23,28�ı߽�����(y����Լ��)
        j = j + 1 ;
        gBC1(j,1) = i ;  % ����
        gBC1(j,2) = 2 ;  % ���ɶȺ�
        gBC1(j,3) = 0 ;  % Լ��ֵ
    end
    j = 6 ;
    for i=5:5:25     %  ���5,10,15,20,25�ı߽�����(y����Լ��)
        j = j + 1;
        gBC1(j,1) = i ;  % ����
        gBC1(j,2) = 2 ;  % ���ɶȺ�
        gBC1(j,3) = 0 ;  % Լ��ֵ
    end
    
    % бԼ������
    bc3_number = 17 ;
    gBC3 = zeros( bc3_number, 4 ) ;
    j = 0 ;
    for i=1:5:26      %  ���1,6,11,16,21,26�ı߽�����(y'����Լ��)
        j = j + 1 ;
        gBC3(j,1) = i ;     % ����
        gBC3(j,2) = 2 ;     % ���ɶȺ�
        gBC3(j,3) = 0 ;     % Լ��ֵ
        gBC3(j,4) = theta ; % ��б�Ƕ�
    end
    j = 6 ;
    for i=4:5:24     %  ���4,9,14,19,24�ı߽�����(y'����Լ��)
        j = j + 1;
        gBC3(j,1) = i ;     % ����
        gBC3(j,2) = 2 ;     % ���ɶȺ�
        gBC3(j,3) = 0 ;     % Լ��ֵ
        gBC3(j,4) = theta ; % ��б�Ƕ�
    end
    j = 11 ;
    for i=2:5:27 
        j = j+1 ;
        gBC3(j,1) = i ;       % ����
        gBC3(j,2) = 2 ;       % ���ɶȺ�
        gBC3(j,3) = 0 ;       % Լ��ֵ
        gBC3(j,4) = theta/2 ; % ��б�Ƕ�
    end
    
    % �ֲ�ѹ��
    df_number = 1 ;
    gDF = [3 1 2 120e6 120e6 120e6] ;
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

    global gNode gElement gMaterial gBC1 gBC3 gDF gK gDelta1 gDelta gElementStress

    % step1. ��������նȾ���ͽڵ�������
    [node_number,dummy] = size( gNode ) ;
    gK = sparse( node_number * 2, node_number * 2 ) ;
    f = sparse( node_number * 2, 2 ) ;    

    % step2. ���㵥Ԫ�նȾ��󣬲����ɵ�����նȾ�����
    [element_number,dummy] = size( gElement ) ;
    for ie=1:1:element_number
        disp( sprintf(  '���㵥Ԫ�նȾ��󲢼��ɣ���ǰ��Ԫ: %d', ie  ) ) ;
        k = StiffnessMatrix( ie ) ;
        AssembleStiffnessMatrix( ie, k ) ;
    end
    
    % step3. ���㵥Ԫ�������ĵ�Ч������������ɵ���������������
    for ie=1:1:element_number
        evf = EquivalentVolumeForce( ie ) ;
        node = gElement( ie, 1:8 ) ;
        f( (node-1)*2+1,2 ) = f( (node-1)*2+1,2 ) + evf(1:2:15) ;
        f( (node-1)*2+2,2 ) = f( (node-1)*2+2,2 ) + evf(2:2:16) ;
    end

    % step4. ����ֲ�ѹ���ĵ�Ч�ڵ����������ɵ�����ڵ���������
    [df_number, dummy] = size( gDF ) ;
    for idf = 1:1:df_number
        enf = EquivalentDistPressure( gDF(idf,1:3), gDF(idf,4:6) ) ;
        i = gDF(idf, 1) ;
        j = gDF(idf, 2) ;
        m = gDF(idf, 3) ;
        f( (i-1)*2+1:(i-1)*2+2, 1 ) = f( (i-1)*2+1:(i-1)*2+2, 1 ) + enf( 1:2 ) ;
        f( (j-1)*2+1:(j-1)*2+2, 1 ) = f( (j-1)*2+1:(j-1)*2+2, 1 ) + enf( 3:4 ) ;
        f( (m-1)*2+1:(m-1)*2+2, 1 ) = f( (m-1)*2+1:(m-1)*2+2, 1 ) + enf( 5:6 ) ;
    end
  
    % step5. ����бԼ���߽�����
    [bc3_number,dummy] = size(gBC3) ;
    for ibc=1:1:bc3_number
        n = gBC3( ibc, 1 ) ;
        theta = gBC3( ibc, 4 ) ;
        c = cos(theta) ;
        s = sin(theta) ;
        T = [ c s
             -s c ] ;
        gK((n-1)*2+1:(n-1)*2+2,:) = T*gK((n-1)*2+1:(n-1)*2+2,:) ;
        gK(:,(n-1)*2+1:(n-1)*2+2) = gK(:,(n-1)*2+1:(n-1)*2+2)*transpose(T) ;
        f((n-1)*2+1:(n-1)*2+2,:) = T*f((n-1)*2+1:(n-1)*2+2,:) ;
        gBC1 = [gBC1; gBC3(ibc,1:3)] ;
    end
    
    % step6. �����һ��Լ���������޸ĸնȾ���ͽڵ������������ó˴�����
    [bc1_number,dummy] = size( gBC1 ) ;
    for ibc=1:1:bc1_number
        n = gBC1(ibc, 1 ) ;
        d = gBC1(ibc, 2 ) ;
        m = (n-1)*2 + d ;
        f(m,:) = gBC1(ibc, 3)* gK(m,m) * 1e15 ;
        gK(m,m) = gK(m,m) * 1e15 ;
    end
    
    % step7. ��ⷽ���飬�õ��ڵ�λ������
    gDelta1 = gK \ f ;
    
    % step8. ����бԼ���Ľ��λ��ת������������
    gDelta = gDelta1 ;
    for ibc=1:bc3_number
        n = gBC3( ibc, 1 ) ;
        theta = gBC3( ibc, 4 ) ;
        c = cos(theta) ;
        s = sin(theta) ;
        T = [ c s
             -s c ] ;
        gDelta( (n-1)*2+1:n*2,:) = transpose(T)*gDelta( (n-1)*2+1:n*2,: ) ;
    end
    
    % step9. ����ÿ����Ԫ�Ľ��Ӧ��
    gElementStress = zeros( element_number, 8, 3, 2 ) ;
    delta = zeros( 16, 2 ) ;
    for ie = 1:element_number
        xi  = [ -1   1   1  -1   0   1   0  -1 ] ;
        eta = [ -1  -1   1   1  -1   0   1   0 ] ;
        for n=1:8
            B = MatrixB( ie, xi(n), eta(n) ) ;
            D = MatrixD( ie, 1 ) ;
            delta(1:2:15, :) = gDelta( gElement(ie,1:8)*2-1, : ) ;
            delta(2:2:16, :) = gDelta( gElement(ie,1:8)*2, : ) ;
            sigma = D*B*delta ;
            gElementStress( ie, n, :, : ) = sigma ;
        end
    end
return

function k = StiffnessMatrix( ie )   
%  ����ƽ��Ӧ��Ȳ�����Ԫ�ĸնȾ���
%  ���������
%     ie -- ��Ԫ��
%  ����ֵ��
%     k  -- ��Ԫ�նȾ���
%  ˵����
%     �ø�˹���ַ����ƽ��Ȳ�����Ԫ�ĸնȾ���

    k = zeros( 16, 16 ) ;
    D = MatrixD( ie, 1 ) ;
    % 3 x 3  ��˹���ֵ��Ȩϵ��
    %x = [-0.774596669241483,              0.0,0.774596669241483] ;   
    %w = [ 0.555555555555556,0.888888888888889,0.555555555555556] ;  
    % 2 x 2  ��˹���ֵ��Ȩϵ��
    x = [-0.577350269189626, 0.577350269189626] ;                    
    w = [1, 1] ;                                                     
    for i=1:1:length(x)
        for j=1:1:length(x)
            B = MatrixB( ie, x(i), x(j) ) ;
            J = Jacobi( ie, x(i), x(j) ) ;
            k = k + w(i)*w(j)*transpose(B)*D*B*det(J) ;   
        end
    end
return

function D = MatrixD( ie, opt )
%  ���㵥Ԫ�ĵ��Ծ���D
%  ���������
%     ie --------- ��Ԫ��
%    opt --------- ����ѡ��
%                   1 ==>  ƽ��Ӧ��
%                   2 ==>  ƽ��Ӧ��
%  ����ֵ��
%     D  --------- ���Ծ���D

    global gElement gMaterial 
    E  = gMaterial( gElement(ie, 9), 1 ) ;   % ����ģ��
    mu = gMaterial( gElement(ie, 9), 2 ) ;   % ���ɱ� 
    if opt == 1   % ƽ��Ӧ���ĵ��Գ���
        A1 = mu ;                                
        A2 = (1-mu)/2 ;                          
        A3  = E/(1-mu^2) ;                      
    else          % ƽ��Ӧ��ĵ��Գ���
       A1 = mu/(1-mu) ;                        
        A2 = (1-2*mu)/2/(1-mu) ;                
        A3 = E*(1-mu)/(1+mu)/(1-2*mu) ;         
    end
    D = A3* [  1  A1   0
              A1   1   0
               0   0  A2] ;
return

function B = MatrixB( ie, xi, eta )
%  ���㵥Ԫ��Ӧ�����B
%  ���������
%     ie --------- ��Ԫ��
%     xi,eta ----- �ֲ�����  
%  ����ֵ��
%     B  --------- �ھֲ����괦��Ӧ�����B

    [N_x,N_y] = N_xy( ie, xi, eta ); 
    B = zeros( 3, 16 ) ;
    for i=1:1:8
        B(1:3,(2*i-1):2*i) = [ N_x(i)        0      
                                    0   N_y(i)
                               N_y(i),  N_x(i)]; 
    end
return

function [N_x, N_y] = N_xy( ie, xi, eta )
%  �����κ�������������ĵ���
%  ���������
%     ie --------- ��Ԫ��
%     xi,eta ----- �ֲ�����  
%  ����ֵ��
%     N_x  ------- �ھֲ����괦���κ�����x����ĵ���
%     N_y  ------- �ھֲ����괦���κ�����y����ĵ���

    J = Jacobi( ie, xi, eta ) ;
    [N_xi,N_eta] = N_xieta( ie, xi, eta ) ;
    A=inv(J)*[N_xi;N_eta] ;
    N_x = A(1,:) ;
    N_y = A(2,:) ;
return

function [N_xi, N_eta] = N_xieta( ie, xi, eta )
%  �����κ����Ծֲ�����ĵ���
%  ���������
%     ie --------- ��Ԫ��
%     xi,eta ----- �ֲ�����  
%  ����ֵ��
%     N_xi   ------- �ھֲ����괦���κ�����xi����ĵ���
%     N_eta  ------- �ھֲ����괦���κ�����eta����ĵ���

    x = [ -1, 1, 1, -1 ] ;
    e = [ -1, -1, 1, 1 ] ;
    N_xi  = zeros( 1, 8 ) ;
    N_eta = zeros( 1, 8 ) ;

    N_xi( 5 )  = xi*(eta-1) ;
    N_eta( 5 ) = 0.5*(xi^2-1) ;
    N_xi( 6 )  = 0.5*(1-eta^2) ;
    N_eta( 6 ) = -eta*(xi+1) ;
    N_xi( 7 )  = -xi*(eta+1) ;
    N_eta( 7 ) = 0.5*(1-xi^2) ;
    N_xi( 8 )  = 0.5*(eta^2-1) ;
    N_eta( 8 ) = eta*(xi-1) ;

    N_xi(1)  = x(1)*(1+e(1)*eta)/4 - 0.5*( N_xi(5)  + N_xi(8) );
    N_eta(1) = e(1)*(1+x(1)*xi)/4  - 0.5*( N_eta(5) + N_eta(8) ) ;
    N_xi(2)  = x(2)*(1+e(2)*eta)/4 - 0.5*( N_xi(5)  + N_xi(6) );
    N_eta(2) = e(2)*(1+x(2)*xi)/4  - 0.5*( N_eta(5) + N_eta(6) ) ;
    N_xi(3)  = x(3)*(1+e(3)*eta)/4 - 0.5*( N_xi(6)  + N_xi(7) );
    N_eta(3) = e(3)*(1+x(3)*xi)/4  - 0.5*( N_eta(6) + N_eta(7) ) ;
    N_xi(4)  = x(4)*(1+e(4)*eta)/4 - 0.5*( N_xi(7)  + N_xi(8) );
    N_eta(4) = e(4)*(1+x(4)*xi)/4  - 0.5*( N_eta(7) + N_eta(8) ) ;
return

function J = Jacobi( ie, xi, eta )
%  �����ſ˱Ⱦ���
%  ���������
%     ie --------- ��Ԫ��
%     xi,eta ----- �ֲ�����  
%  ����ֵ��
%     J   ------- �ھֲ�����(xi,eta)�����ſ˱Ⱦ���
    global gNode gElement
    x = gNode(gElement(ie,1:8),1) ;
    y = gNode(gElement(ie,1:8),2) ;
    [N_xi,N_eta] = N_xieta( ie, xi, eta ) ;
    x_xi  = N_xi  * x ;
    x_eta = N_eta * x ;
    y_xi  = N_xi  * y ;
    y_eta = N_eta * y ;
    J = [ x_xi, y_xi; x_eta, y_eta ];
return

function AssembleStiffnessMatrix( ie, k )
%  �ѵ�Ԫ�նȾ��󼯳ɵ�����նȾ���
%  �������:
%      ie  --- ��Ԫ��
%      k   --- ��Ԫ�նȾ���
%  ����ֵ:
%      ��
    global gElement gK
    for i=1:1:8
        for j=1:1:8
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

function edf = EquivalentDistPressure( node, pressure )
%   ����ֲ�ѹ���ĵ�Ч�ڵ���
%   �������:
%      node ----------  ����
%      pressure ------  �����Ŷ�Ӧ��ѹ��ֵ
%   ����ֵ:
%      edf -----------  ��Ч�ڵ�������
    global gNode
    x = gNode( node, 1 ) ;
    y = gNode( node, 2 ) ;
    xi = x(1) ; xj = x(2) ; xk = x(3) ;
    yi = y(1) ; yj = y(2) ; yk = y(3) ;
    sigi = pressure(1) ;
    sigj = pressure(2) ;
    sigk = pressure(3) ;
    
    X1 = -10*xi-2*xj+12*xk ;
    Y1 = -10*yi-2*yj+12*yk ;
    X2 = xi - xj ;
    Y2 = yi - yj ;
    X3 = -6*xi-2*xj+8*xk ;
    Y3 = -6*yi-2*yj+8*yk ;
    X4 = 2*xi+10*xj-12*xk ;
    Y4 = 2*yi+10*yj-12*yk ;
    X5 = 2*xi+6*xj-8*xk ;
    Y5 = 2*yi+6*yj-8*yk ;
    X6 = -16*xi+16*xj ;
    Y6 = -16*yi+16*yj ;
    
    XY = [  Y1  Y2  Y3
           -X1 -X2 -X3
            Y2  Y4  Y5
           -X2 -X4 -X5
            Y3  Y5  Y6
           -X3 -X5 -X6 ] ;
    sig = [sigi; sigj; sigk ] ;
    edf = 1/30*XY*sig ;
return

function evf = EquivalentVolumeForce( ie )
%   �����������ĵ�Ч�ڵ���
%   �������:
%      ie ----------  ��Ԫ��
%   ����ֵ:
%      evf ---------  ��Ч�ڵ�������
    global gNode gElement gMaterial
    
    evf = zeros( 16, 1 ) ;
    omega = 2094.0 ;   % ��ת���ٶ�
    ro = gMaterial( gElement( ie, 9 ), 3 ) ;
    x = gNode( gElement(ie,1:8), 1 ) ;
    y = gNode( gElement(ie,1:8), 2 ) ;
    xi = [-0.774596669241483,0.0,              0.774596669241483] ;
    w  = [ 0.555555555555556,0.888888888888889,0.555555555555556] ;
    for i=1:length(xi)
        for j=1:length(xi)
            J = Jacobi( ie, xi(i), xi(j) ) ;
            detJ = det(J); 
            N = ShapeFunction( xi(i), xi(j) ) ;
            evf(1:2:15) = evf(1:2:15) + N'*N*x*detJ*w(i)*w(j) ;
            evf(2:2:16) = evf(2:2:16) + N'*N*y*detJ*w(i)*w(j) ;
        end
    end
    evf = evf * ro * omega^2 ;
return

function N = ShapeFunction( xi, eta )
%   �����κ�����ֵ
%   �������:
%      ie ----------  ��Ԫ��
%      xi, eta -----  ��Ԫ�ھֲ�����
%   ����ֵ:
%      N -----------  �κ�����ֵ

    N5 = ( eta - 1 ) * ( xi^2 - 1 ) / 2 ;
    N6 = ( xi + 1 ) * ( 1 - eta^2 ) / 2 ;
    N7 = ( eta + 1 ) * ( 1 - xi^2 ) / 2 ;
    N8 = ( xi - 1 ) * ( eta^2 - 1 ) / 2 ;
    N1 = ( 1 - xi ) * ( 1 - eta ) / 4 - 0.5 * ( N8 + N5 ) ;
    N2 = ( 1 + xi ) * ( 1 - eta ) / 4 - 0.5 * ( N5 + N6 ) ;
    N3 = ( 1 + xi ) * ( 1 + eta ) / 4 - 0.5 * ( N6 + N7 ) ;
    N4 = ( 1 - xi ) * ( 1 + eta ) / 4 - 0.5 * ( N7 + N8 ) ;
    N = [ N1 N2 N3 N4 N5 N6 N7 N8 ]; 
return

function SaveResults( file_out ) 
%  ���������
%  ���������
%     ��
%  ����ֵ��
%     ��
    global gNode gElement gMaterial gBC1 gBC3 gDF gK gDelta gDelta1 gElementStress
    save( file_out, 'gNode', 'gElement', 'gMaterial', 'gBC1', 'gBC3', ...
          'gDF', 'gDelta', 'gDelta1', 'gElementStress' ) ;
return

function DisplayResults
%  ��ʾ��������������������Ƚ�
%  ���������
%     ��
%  ����ֵ��
%     ��

%  ˵��:
%    1������ѹ����Գ�ԲͲƽ��Ӧ������Ľ�������(�ο�����[1])
%           ����λ��:  ur = -C1/E*(1+mu)/r + 2*C2/E*(1-mu)*r
%                      C1 = -a^2*b^2*p/(b^2-a^2)
%                      C2 = p*a^2/2/(b^2-a^2)
%           ����Ӧ��:  sr = a^2*p/(b^2-a^2)*(1-b^2/r^2)
%           ����Ӧ��:  st = a^2*p/(b^2-a^2)*(1+b^2/r^2)
%    2�����Ƚ��ٶ�w��ת��ԲͲƽ��Ӧ������Ľ�����(�ο�����[2])
%           ����λ��:  ur = 1/E*((1-mu)*C3*r - (1+mu)*C4/r - (1-mu^2)/8*ro*w^2*r^3 ) 
%                      C3 = (3+mu)/8*ro*w^2*(b^2+a^2)
%                      C4 = -(3+mu)/8*ro*w^2*a^2*b^2
%           ����Ӧ��:  sr = (3+mu)/8*ro*w^2*(b^2+a^2-a^2*b^2/r^2-r^2)
%           ����Ӧ��:  st = (3+mu)/8*ro*w^2*(b^2+a^2+a^2*b^2/r^2-(1+3*mu)/(3+mu)*r^2)
%     
%   ����˵��
%        E  ----- ����ģ��
%       mu  ----- ���ɱ�
%        a  ----- �ڰ뾶
%        b  ----- ��뾶
%        r  ----- ��������
%        p  ----- ��ѹ
%       ro  ----- �ܶ�
%        w  ----- ��ת���ٶ�(rad/s)
%  
%   �ο�����
%     [1] Ǯΰ�� Ҷ��Ԫ, 1956, ������ѧ, pp.230-231
%     [2] Timoshenko S. P., Goodier J. N., 1970 Theory of elasticity. New York: McGraw-Hill, pp.80-83

    global gNode gElement gMaterial gDelta gDelta1 gElementStress

    E  = gMaterial( 1, 1 ) ;
    mu = gMaterial( 1, 2 ) ;
    ro = gMaterial( 1, 3 ) ;
    a = 0.05 ;
    b = 0.10 ;
    p = 120e6 ;
    w = 2094 ;
    r = 0.05:0.005:0.1 ;

    %  ��������ѹԲͲ��λ�ƺ�Ӧ���Ľ�����
    C1 = -a^2*b^2*p/(b^2-a^2) ;
    C2 = p*a^2/2/(b^2-a^2) ;
    ur1 = -C1/E*(1+mu)./r + 2*C2/E*(1-mu)*r ;
    sr1 = a^2*p/(b^2-a^2)*(1-b^2./r.^2) ;
    st1 = a^2*p/(b^2-a^2)*(1+b^2./r.^2) ;
    
    %  ������תԲ�̵�λ�ƺ�Ӧ���Ľ�����
    C3 = (3+mu)/8*ro*w^2*(b^2+a^2) ;
    C4 = -(3+mu)/8*ro*w^2*a^2*b^2 ;
    ur2 = 1/E*((1-mu)*C3*r - (1+mu)*C4./r - (1-mu^2)/8*ro*w^2*r.^3 ) ;
    sr2 = (3+mu)/8*ro*w^2*(b^2+a^2-a^2*b^2./r.^2-r.^2) ;
    st2 = (3+mu)/8*ro*w^2*(b^2+a^2+a^2*b^2./r.^2-(1+3*mu)/(3+mu)*r.^2) ;
    
    %  ��ȡ����Ԫ�ļ�����
    node = [3, 5, 8, 10, 13, 15, 18, 20, 23, 25, 28]; 
    node_num = zeros( size(node) ) ;
    ur1_fem = gDelta( node*2-1, 1 ) ;
    ur2_fem = gDelta( node*2-1, 2 ) ;
    [element_number,dummy] = size( gElement ) ;
    sr1_fem = zeros( size( sr1 ) ) ;
    st1_fem = zeros( size( st1 ) ) ;
    sr2_fem = zeros( size( sr2 ) ) ;
    st2_fem = zeros( size( st2 ) ) ;
    for n=1:length(node)
        for ie=1:element_number
            index = find( gElement( ie, 1:8 ) == node( n ) ) ;
            if ~isempty(index)    
                sr1_fem( n ) = sr1_fem( n ) + gElementStress( ie, index, 1, 1 ) ;
                st1_fem( n ) = st1_fem( n ) + gElementStress( ie, index, 2, 1 ) ;
                sr2_fem( n ) = sr2_fem( n ) + gElementStress( ie, index, 1, 2 ) ;
                st2_fem( n ) = st2_fem( n ) + gElementStress( ie, index, 2, 2 ) ;
                node_num(n) = node_num(n) + 1 ;
            end
        end
    end
    sr1_fem = sr1_fem ./ node_num ;
    st1_fem = st1_fem ./ node_num ;
    sr2_fem = sr2_fem ./ node_num ;
    st2_fem = st2_fem ./ node_num ;
    
    fprintf( '\n\n\n' ) ;
    fprintf( '                        ��1 ����λ�� (mm)\n' ) ;
    fprintf( '==================================================================\n' ) ;
    fprintf( '����  ����         ������ѹ                     ����������\n' ) ;
    fprintf( '        (mm)    ����Ԫ��      ������        ����Ԫ��     ������\n' ) ;
    fprintf( '------------------------------------------------------------------\n' ) ;
    for i=1:length(node)
        fprintf( '%4d   %4.0f     %.6f     %.6f       %.6f     %.6f\n', ...
             node(i), r(i)*1000, ur1_fem(i)*1000, ur1(i)*1000, ur2_fem(i)*1000, ur2(i)*1000 ) ;
    end
    fprintf( '==================================================================\n\n\n' ) ;
    

    fprintf( '                      ��2 ����Ӧ�� (MPa)\n' ) ;
    fprintf( '==================================================================\n' ) ;
    fprintf( '����  ����         ������ѹ                     ����������\n' ) ;
    fprintf( '        (mm)    ����Ԫ��      ������        ����Ԫ��     ������\n' ) ;
    fprintf( '------------------------------------------------------------------\n' ) ;
    for i=1:length(node)
        fprintf( '%4d   %4.0f    %8.2f    %8.2f       %8.2f     %8.2f\n', ...
             node(i), r(i)*1000, sr1_fem(i)/1e6, sr1(i)/1e6, sr2_fem(i)/1e6, sr2(i)/1e6 ) ;
    end
    fprintf( '==================================================================\n\n\n' ) ;

    
    fprintf( '                      ��3 ����Ӧ�� (MPa)\n' ) ;
    fprintf( '==================================================================\n' ) ;
    fprintf( '����  ����         ������ѹ                     ����������\n' ) ;
    fprintf( '        (mm)    ����Ԫ��      ������        ����Ԫ��     ������\n' ) ;
    fprintf( '------------------------------------------------------------------\n' ) ;
    for i=1:length(node)
        fprintf( '%4d   %4.0f    %8.2f    %8.2f       %8.2f     %8.2f\n', ...
             node(i), r(i)*1000, st1_fem(i)/1e6, st1(i)/1e6, st2_fem(i)/1e6, st2(i)/1e6 ) ;
    end
    fprintf( '==================================================================\n\n\n' ) ;
return