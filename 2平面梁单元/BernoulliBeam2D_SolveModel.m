function BernoulliBeam2D_SolveModel
%      �ú����������Ԫģ�ͣ���������
%        1. ���㵥Ԫ�նȾ��󣬼�������նȾ���
%        2. ���㵥Ԫ�ĵ�Ч�ڵ�������������ڵ�������
%        3. ����Լ���������޸�����նȾ���ͽڵ�������
%        4. ��ⷽ���飬�õ�����ڵ�λ������

    global gNode gElement gMaterial gBC1 gNF gDF gK gDelta

    % step1. ��������նȾ���ͽڵ�������
    [node_number,dummy] = size( gNode ) ;
%     sparse ����ϡ�����,gDelta = gK \ fҲ��ϡ�����
    gK = sparse( node_number * 3, node_number * 3 ) ;
    f = sparse( node_number * 3, 1 ) ;

    % step2. ���㵥Ԫ�նȾ��󣬲����ɵ�����նȾ�����
    [element_number,dummy] = size( gElement ) ;
    for ie=1:1:element_number
        k = BernoulliBeam2D_Stiffness( ie, 1 ) ;
        BernoulliBeam2D_Assembly( ie, k ) ;
    end

    % step3. �Ѽ�����ֱ�Ӽ��ɵ�����ڵ���������
    [nf_number, dummy] = size( gNF ) ;
    for inf=1:1:nf_number
        n = gNF( inf, 1 ) ;
        d = gNF( inf, 2 ) ;
        f( (n-1)*3 + d ) = gNF( inf, 3 ) ;
    end

    % step4. ����ֲ����ĵ�Ч�ڵ��������ɵ�����ڵ���������
    [df_number, dummy] = size( gDF ) ;
    for idf = 1:1:df_number
        enf = BernoulliBeam2D_EquivalentNodeForce( gDF(idf,1), gDF(idf, 2), gDF( idf, 3), gDF( idf, 4 ) ) ;
        i = gElement( gDF(idf,1), 1 ) ;
        j = gElement( gDF(idf,1), 2 ) ;
        f( (i-1)*3+1 : (i-1)*3+3 ) = f( (i-1)*3+1 : (i-1)*3+3 ) + enf( 1:3 ) ;
        f( (j-1)*3+1 : (j-1)*3+3 ) = f( (j-1)*3+1 : (j-1)*3+3 ) + enf( 4:6 ) ;
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
return