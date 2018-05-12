function BernoulliBeam2D_SolveModel
%      �ú����������Ԫģ�ͣ���������
%        1. ���㵥Ԫ�նȾ��󣬼�������նȾ���
%        2. ���㵥Ԫ�ĵ�Ч�ڵ�������������ڵ�������
%        3. ����Լ���������޸�����նȾ���ͽڵ�������
%        4. ��ⷽ���飬�õ�����ڵ�λ������

    global gNode gElement gMaterial gBC1 gNF gDF gK gDelta

    % step1. ��������նȾ���ͽڵ�������
    [node_number,~] = size( gNode ) ;
%     sparse ����ϡ�����,gDelta = gK \ fҲ��ϡ�����
    gK = sparse( node_number * 3, node_number * 3 ) ;
    f = sparse( node_number * 3, 1 ) ;

    % step2. ���㵥Ԫ�նȾ��󣬲����ɵ�����նȾ�����
    [element_number,~] = size( gElement ) ;
    for ie=1:1:element_number
        k = BernoulliBeam2D_Stiffness( ie, 1 ) ;
        BernoulliBeam2D_Assembly( ie, k ) ;
    end

    % step3. �Ѽ�����ֱ�Ӽ��ɵ�����ڵ���������
    [nf_number, ~] = size( gNF ) ;
    for inf=1:1:nf_number
        n = gNF( inf, 1 ) ;
        d = gNF( inf, 2 ) ;
        f( (n-1)*3 + d ) = gNF( inf, 3 ) ;
    end

    % step4. ����ֲ����ĵ�Ч�ڵ��������ɵ�����ڵ���������
    [df_number, ~] = size( gDF ) ;
    for idf = 1:1:df_number
        enf = BernoulliBeam2D_EquivalentNodeForce( gDF(idf,1), gDF(idf, 2), gDF( idf, 3), gDF( idf, 4 ) ) ;
        i = gElement( gDF(idf,1), 1 ) ;
        j = gElement( gDF(idf,1), 2 ) ;
        f( (i-1)*3+1 : (i-1)*3+3 ) = f( (i-1)*3+1 : (i-1)*3+3 ) + enf( 1:3 ) ;
        f( (j-1)*3+1 : (j-1)*3+3 ) = f( (j-1)*3+1 : (j-1)*3+3 ) + enf( 4:6 ) ;
    end
  
    % step5. ����Լ���������޸ĸնȾ���ͽڵ������������ó˴�����
    [bc_number,~] = size( gBC1 ) ;
    for ibc=1:1:bc_number
        n = gBC1(ibc, 1 ) ;
        d = gBC1(ibc, 2 ) ;
        m = (n-1)*3 + d ;
        f(m) = gBC1(ibc, 3)* gK(m,m) * 1e15 ;
        gK(m,m) = gK(m,m) * 1e15 ;
    end

    % step 6. ��ⷽ���飬�õ��ڵ�λ������
    gDelta = gK \ f ;
return

function k=BernoulliBeam2D_Stiffness(ie,icoord)
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
        T = BernoulliBeam2D_TransformMatrix( ie ) ;
        k = T*k*transpose(T) ;
    end
return

function BernoulliBeam2D_Assembly( ie, k )
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

function T = BernoulliBeam2D_TransformMatrix( ie )
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

function enf = BernoulliBeam2D_EquivalentNodeForce( ie, p1, p2, idof )
%   �������Էֲ����صĵ�Ч�ڵ���
%   �������:
%      ie  -----  ��Ԫ��
%      p1  -----  ��һ���ڵ��ϵķֲ�������ֵ
%      p2  -----  �ڶ����ڵ��ϵķֲ�������ֵ
%      idof  ---  �ֲ��������࣬�����������漸��
%                  1 ---  �ֲ�������
%                  2 ---  �ֲ�������
%                  3 ---  �ֲ����
%   ����ֵ:
%      enf -----  ��������ϵ�µ�Ч�ڵ�������
    global gElement gNode
    enf = zeros( 6, 1 ) ; % ���� 6x1 �ĵ�Ч�ڵ�������
    xi = gNode( gElement( ie, 1 ), 1 ) ;
    yi = gNode( gElement( ie, 1 ), 2 ) ;
    xj = gNode( gElement( ie, 2 ), 1 ) ;
    yj = gNode( gElement( ie, 2 ), 2 ) ;
    L = sqrt( (xj-xi)^2 + (yj-yi)^2 ) ;
    switch idof 
    case 1     %  �ֲ������� 
        enf( 1 ) = (2*p1+p2)*L/6 ;
        enf( 4 ) = (p1+2*p2)*L/6 ;
    case 2     %  �ֲ�������
        enf( 2 ) = (7*p1+3*p2)*L/20 ;
        enf( 3 ) = (3*p1+2*p2)*L^2/60 ;
        enf( 5 ) = (3*p1+7*p2)*L/20 ;
        enf( 6 ) = -(2*p1+3*p2)*L^2/60 ;
    case 3     %  �ֲ����
        enf( 2 ) = -(p1+p2)/2 ;
        enf( 3 ) = (p1-p2)*L/12 ;
        enf( 5 ) = (p1+p2)/2 ;
        enf( 6 ) = -(p1-p2)*L/12 ;
    otherwise
        disp( sprintf( '�ֲ�����������󣬵�Ԫ��:%d',ie ) ) ;
    end
    
    T = BernoulliBeam2D_TransformMatrix( ie );%  ���㵥Ԫ��ת������
    enf = T * enf;%  �ѵ�Ч�ڵ���ת��������������
return