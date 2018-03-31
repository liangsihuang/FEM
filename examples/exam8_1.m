function exam8_1
% ������Ϊ�ڰ��µĵ�һ������������ƽ������Ԫ�������������߹�������������
%      ��������� ��
%      �������� ǰ3����Ƶ�ʼ�����Ӧ������ 

% ����ȫ�ֱ���
%      gNode ------ �ڵ�����
%      gElement --- ��Ԫ����
%      gMaterial -- ��������
%      gBC1 ------- ��һ��Լ������
%      gK --------- ����նȾ���
%      gDelta ----- ����ڵ�����

    PlaneFrameModel ;      % ��������Ԫģ��
    SolveModel ;           % �������Ԫģ��
    DisplayResults ;       % ��ʾ������
return ;

function PlaneFrameModel
%  ����ƽ���ϵ������Ԫģ��
%  ���������
%      ��
%  ����ֵ��
%      ��
%  ˵����
%      �ú�������ƽ���ϵ������Ԫģ�����ݣ�
%        gNode ------- �ڵ㶨��
%        gElement ---- ��Ԫ����
%        gMaterial --- ���϶��壬��������ģ�������Ľ���������Ŀ�����Ծ�
%        gBC --------- Լ������

    global gNode gElement gMaterial gBC1

    % ���������߹��ļ�������
    L = 60 ;               %  ����羶     
    f = 7.5 ;              %  ����ʸ��
    
    n = 100 ;              %  ��Ԫ��Ŀ
    x = -L/2:L/n:L/2 ;     %  ����x����
    a = f/L^2*4 ;
    y = - a * x.^2 ;       %  ����y����

    % �ڵ�����
    gNode = [x'  y'] ;
    
    % ��Ԫ����
    gElement = zeros( n, 3 ) ;
    for i=1:n
        gElement( i, : ) = [ i, i+1, 1 ] ;
    end
    
    % �������� 
    %           ����ģ��   ������Ծ�   �����   �ܶ�
    gMaterial = [2.06e11,  0.03622,   0.0815,  1435.2/0.0815];   %  ���� 1

    % ��һ��Լ������
    %     �ڵ��   ���ɶȺ�    Լ��ֵ
    gBC1 = [ 1,        1,        0.0
             1,        2,        0.0
             n+1,      1,        0.0
             n+1,      2,        0.0] ;
return

function SolveModel
%  �������Ԫģ��
%  ���������
%     ��
%  ����ֵ��
%     ��
%  ˵����
%      �ú����������Ԫģ�ͣ���������
%        1. ���㵥Ԫ�ĸնȺ��������󣬼�������նȺ���������
%        2. ����Լ���������޸�����նȾ���
%        3. �������ֵ����

    global gNode gElement gMaterial gBC1 gK gM gEigValue gEigVector

    % step1. ��������նȾ���ͽڵ�������
    [node_number,dummy] = size( gNode ) ;
    gK = sparse( node_number * 3, node_number * 3 ) ;
    gM = sparse( node_number * 3, node_number * 3 ) ;

    % step2. ���㵥Ԫ�նȺ��������󣬲����ɵ�����նȺ�����������
    [element_number,dummy] = size( gElement ) ;
    for ie=1:1:element_number
        k = StiffnessMatrix( ie ) ;
        m = MassMatrix( ie ) ; 
        AssembleGlobalMatrix( ie, k, m ) ;
    end

    % step3. �����һ��Լ���������޸ĸնȾ�����������󡣣����û��л��з���
    [bc1_number,dummy] = size( gBC1 ) ;
    w2max = max( diag(gK)./diag(gM) ) ;
    for ibc=1:1:bc1_number
        n = gBC1(ibc, 1 ) ;
        d = gBC1(ibc, 2 ) ;
        m = (n-1)*3 + d ;
        gK(:,m) = zeros( node_number*3, 1 ) ;
        gK(m,:) = zeros( 1, node_number*3 ) ;
        gK(m,m) = 1;
        gM(:,m) = zeros( node_number*3, 1 ) ;
        gM(m,:) = zeros( 1, node_number*3 ) ;
        gM(m,m) = gK(m,m)/w2max/1e10 ;
    end
    
    % step4. �������ֵ����
    % step4.1Ϊ��ʹ�նȾ������������Գƣ��ڼ���ʱ��������������
    for i=1:node_number*3
        for j=i:node_number*3
            gK(j,i) = gK(i,j) ;
            gM(j,i) = gM(i,j) ;
        end
    end
    
    % step4.2 ����ǰ6������ֵ����������
    [gEigVector, gEigValue] = eigs(gK, gM, 3, 'SM' ) ; 
    
    % step4.3 �޸�������������Լ�������ɶ�
    for ibc=1:1:bc1_number
        n = gBC1(ibc, 1 ) ;
        d = gBC1(ibc, 2 ) ;
        m = (n-1)*3 + d ;
        gEigVector(m,:) = gBC1(ibc,3) ;
    end
return

function k = StiffnessMatrix( ie )
%  ���㵥Ԫ�նȾ���
%  �������:
%     ie -------  ��Ԫ��
%  ����ֵ:
%     k  ----  ��������ϵ�µĸնȾ���
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
    T = TransformMatrix( ie ) ;
    k = T*k*transpose(T) ;
return

function m = MassMatrix( ie )
%  ���㵥Ԫ��������
%  �������:
%     ie -------  ��Ԫ��
%  ����ֵ:
%     m  ----  ��������ϵ�µ���������
    global gNode gElement gMaterial
    m = zeros( 6, 6 ) ;
    E = gMaterial( gElement(ie, 3), 1 ) ;
    A = gMaterial( gElement(ie, 3), 3 ) ;
    ro = gMaterial( gElement(ie, 3 ), 4 ) ;
    xi = gNode( gElement( ie, 1 ), 1 ) ;
    yi = gNode( gElement( ie, 1 ), 2 ) ;
    xj = gNode( gElement( ie, 2 ), 1 ) ;
    yj = gNode( gElement( ie, 2 ), 2 ) ;
    L = ( (xj-xi)^2 + (yj-yi)^2 )^(1/2) ;
    m = ro*A*L/420*[140      0      0   70      0      0
                      0    156   22*L    0     54   -13*L
                      0   22*L  4*L^2    0   13*L  -3*L^2
                     70      0      0  140      0       0 
                      0     54   13*L    0    156   -22*L
                      0  -13*L -3*L^2    0  -22*L  4*L^2 ] ;
    T = TransformMatrix( ie ) ;
    m = T*m*transpose(T) ;
return

function AssembleGlobalMatrix( ie, ke, me )
%  �ѵ�Ԫ�նȺ��������󼯳ɵ�����նȾ���
%  �������:
%      ie  --- ��Ԫ��
%      ke  --- ��Ԫ�նȾ���
%      me  --- ��Ԫ��������
%  ����ֵ:
%      ��
    global gElement gK gM
    for i=1:1:2
        for j=1:1:2
            for p=1:1:3
                for q =1:1:3
                    m = (i-1)*3+p ;
                    n = (j-1)*3+q ;
                    M = (gElement(ie,i)-1)*3+p ;
                    N = (gElement(ie,j)-1)*3+q ;
                    gK(M,N) = gK(M,N) + ke(m,n) ;
                    gM(M,N) = gM(M,N) + me(m,n) ;
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

function DisplayResults
%  ��ʾ������
%  ���������
%     ��
%  ����ֵ��
%     ��

    global gNode gElement gMaterial gBC1 gEigValue gEigVector

    fre_number = length(diag(gEigValue)) ;
    
    % ��ӡ��������(����)
    fprintf( '\n\n ��һ   ��������(����)  \n' ) ;
    for i=1:fre_number
        fprintf( '----------------') ;
    end
    fprintf( '\n' ) ;
    for i=1:fre_number
        fprintf( '  %6d        ', i ) ;
    end
    fprintf( '\n' ) ;
    for i=1:fre_number
        fprintf( '----------------') ;
    end
    fprintf( '\n' ) ;
    [dof,dummy]=size(gEigVector) ;
    for i=1:dof
        for j=fre_number:-1:1
            fprintf( '%15.7e ', gEigVector(i,j) ) ;
        end
        fprintf( '\n' ) ;
    end
    for i=1:fre_number
        fprintf( '----------------') ;
    end
    fprintf( '\n' ) ;

    % ��ӡ����ֵ
    fprintf( '\n\n\n\n ���   ����ֵ(Ƶ��)�б�  \n' ) ;
    fprintf( '----------------------------------------------------------\n') ;
    fprintf( '   ����        ����ֵ          Ƶ��(Hz)        ԲƵ��(Hz)  \n' ) ;
    fprintf( '----------------------------------------------------------\n') ;
    for i=fre_number:-1:1
        fprintf( '%6d   %15.7e   %15.7e   %15.7e\n', fre_number-i+1, ...
            gEigValue(i,i), sqrt(gEigValue(i,i))/2/pi, sqrt(gEigValue(i,i)) ) ;
    end
    fprintf( '----------------------------------------------------------\n') ;
    
    % ��������ͼ
    for j=fre_number:-1:1
        figure ;
        x = gNode(:,1) ;
        y = gNode(:,2) ;
        dx = gEigVector(1:3:length(x)*3, j ) ;
        dy = gEigVector(2:3:length(x)*3, j ) ;
        factor = max( [max(abs(x))/max(abs(dx)), max(abs(y))/max(abs(dy))] )* 0.05; 
        plot(x,y,'-', x+factor*dx, y+factor*dy, ':') ;
        title( sprintf( '��%d��Ƶ��: %.3f Hz', fre_number-j+1, sqrt(gEigValue(j,j))/2/pi ) ) ;
        axis equal;
        axis off ;
    end
return
