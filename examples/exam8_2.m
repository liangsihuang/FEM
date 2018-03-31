function exam8_2
% ������Ϊ�ڰ��µĵڶ�������������ƽ������Ԫ�������������߹����ڳ�ʼ������
%  �����񶯣�����ʱ�����߽������FFT�任����õ�Ƶ�ʿ���exam8_1.m�Ľ������
%  �Ƚϣ�����֤������Ŀɿ���
%      ��������� ��
%      �������� λ�Ƶ�ʱ�����߼���Ƶ������ͼ 

    PlaneFrameModel ;             % ��������Ԫģ��
    SolveModel ;                  % �������Ԫģ��
    SaveResults('exam8_2.mat') ;  % ���������
    DisplayResults ;              % ��ʾ������
return ;

function PlaneFrameModel
%  ����ƽ���ϵ������Ԫģ��
%  ���������
%      ��
%  ����ֵ��
%      ��
%  ˵����
%      �ú�������ƽ���ϵ������Ԫģ�����ݣ�
%        gNode -------- �ڵ㶨��
%        gElement ----- ��Ԫ����
%        gMaterial ---- ���϶��壬��������ģ�������Ľ���������Ŀ�����Ծ�
%        gBC1 --------- Լ������
%        gDeltaT ------ ʱ�䲽��
%        gTimeEnd ----- �������ʱ��
%        gDisp -------- λ��ʱ����Ӧ
%        gVelo -------- �ٶ�ʱ����Ӧ
%        gAcce -------- ���ٶ�ʱ����Ӧ

    global gNode gElement gMaterial gBC1 gDeltaT gTimeEnd gDisp gVelo gAcce

    % ���������߹��ļ�������
    L = 60 ;               %  ����羶(m)     
    f = 7.5 ;              %  ����ʸ��(m)
    
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
     
    gDeltaT = 0.01 ;
    gTimeEnd = 4096*gDeltaT  ;    % ����ʱ��Ϊ�غ�ͨ������ʱ�������
    timestep = floor(gTimeEnd/gDeltaT) ;

    % ����λ�ƣ��ٶȺͼ��ٶ�
    gDisp = zeros( (n+1)*3, timestep ) ;
    gVelo = zeros( (n+1)*3, timestep ) ;
    gAcce = zeros( (n+1)*3, timestep ) ;
    
    % ��ʼ����
    gDisp(:,1) = zeros( (n+1)*3, 1 ) ;
    gVelo(:,1) = ones( (n+1)*3, 1 ) ;
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
%        2. ��Newmark������ʱ����Ӧ

    global gNode gElement gMaterial gBC1 gK gM gDeltaT gTimeEnd gDisp gVelo gAcce

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

    % step3. ����ʱ����Ӧ(Newmark��)
    % step3.1 ��ʼ����
    gama = 0.5 ;
    beta = 0.25 ;
    C = zeros( size( gK ) ) ;
    [N,N] = size( gK ) ;
    alpha0 = 1/beta/gDeltaT^2 ;
    alpha1 = gama/beta/gDeltaT ;
    alpha2 = 1/beta/gDeltaT ;
    alpha3 = 1/2/beta - 1 ;
    alpha4 = gama/beta - 1 ;
    alpha5 = gDeltaT/2*(gama/beta-2) ;
    alpha6 = gDeltaT*(1-gama) ;
    alpha7 = gama*gDeltaT ;
    K1 = gK + alpha0*gM + alpha1*C;
    timestep = floor(gTimeEnd/gDeltaT) ;
    
    % step3.2 ��K1���б߽���������
    [bc1_number,dummy] = size( gBC1 ) ;
    K1im = zeros(N,bc1_number) ;
    for ibc=1:1:bc1_number
        n = gBC1(ibc, 1 ) ;
        d = gBC1(ibc, 2 ) ;
        m = (n-1)*3 + d ;
        K1im(:,ibc) = K1(:,m) ; 
        K1(:,m) = zeros( node_number*3, 1 ) ;
        K1(m,:) = zeros( 1, node_number*3 ) ;
        K1(m,m) = 1.0 ;
    end
    [KL,KU] = lu(K1) ;   % �������Ƿֽ⣬��ʡ��������ʱ��
    
    % step3.3 �����ʼ���ٶ�
    gAcce(:,1) = gM\(-gK*gDisp(:,1)-C*gVelo(:,1)) ;
    
    % step3.4 ��ÿһ��ʱ�䲽����
    for i=2:1:timestep
        if mod(i,100) == 0
            fprintf( '��ǰʱ�䲽��%d\n', i ) ;
        end
        f1 = gM*(alpha0*gDisp(:,i-1)+alpha2*gVelo(:,i-1)+alpha3*gAcce(:,i-1)) ...
                  + C*(alpha1*gDisp(:,i-1)+alpha4*gVelo(:,i-1)+alpha5*gAcce(:,i-1)) ;
        % ��f1���б߽���������
        [bc1_number,dummy] = size( gBC1 ) ;
        for ibc=1:1:bc1_number
            n = gBC1(ibc, 1 ) ;
            d = gBC1(ibc, 2 ) ;
            m = (n-1)*3 + d ;
            f1 = f1 - gBC1(ibc,3) * K1im(:,ibc) ;
            f1(m) = gBC1(ibc,3) ;
        end
        y = KL\f1 ;
        gDisp(:,i) = KU\y ;
        gAcce(:,i) = alpha0*(gDisp(:,i)-gDisp(:,i-1)) - alpha2*gVelo(:,i-1) - alpha3*gAcce(:,i-1) ;
        gVelo(:,i) = gVelo(:,i-1) + alpha6*gAcce(:,i-1) + alpha7*gAcce(:,i) ;
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

function SaveResults( file_out )
%  ���������
%  ���������
%     ��
%  ����ֵ��
%     ��
    global gNode gElement gMaterial gBC1 gDeltaT gTimeEnd gLoad gLoadVelo gDisp gVelo gAcce
    save( file_out, 'gNode', 'gElement', 'gMaterial', 'gBC1',  ...
          'gDeltaT', 'gTimeEnd', 'gLoad', 'gLoadVelo', 'gDisp', 'gVelo', 'gAcce' ) ;
return

function DisplayResults
%  ��ʾ������
%  ���������
%     ��
%  ����ֵ��
%     ��

    global gNode gElement gMaterial gBC1 gDisp gVelo gAcce gDeltaT gTimeEnd

    % ����ʱ������
    [node_number,dummy] = size(gNode) ;
    t = 0:gDeltaT:gTimeEnd-gDeltaT ;
    d = gDisp((floor(node_number/4)*3)+2,:) ;
    subplot(2,1,1) ;
    plot( t, d ) ;
    title( 'L/4���Ӷ�ʱ������' ) ;
    xlabel( 'ʱ��(s)') ;
    ylabel( '�Ӷ�(m)' ) ;
    
    % ��ʱ�����߽���FFT�任����ȡƵ������
    fd = fft( d ) ;
    df = 1/gTimeEnd ;
    f = (0:length(d)-1)*df ;
    subplot(2,1,2);
    plot(f,abs(fd)) ;
    set(gca,'xlim',[2,10]) ;
    title( 'L/4���Ӷȵ�Ƶ��ͼ' ) ;
    xlabel( 'Ƶ��(Hz)') ;
    ylabel( '��ֵ' ) ;
    
    % ��עƵ�ʷ�ֵ
    fifi1 = diff(abs(fd));
    n = length(fifi1) ;
    d1 = fifi1(1:n-1);
    d2 = fifi1(2:n) ;
    indmax = find( d1.*d2<0 & d1>0 )+1;
    for i=1:length(indmax)
        if f(indmax(i)) > 10 
            break ;
        end
        text( f(indmax(i)+2), abs(fd(indmax(i)))*0.9, sprintf('f=%.3f',f(indmax(i))));
    end
return
