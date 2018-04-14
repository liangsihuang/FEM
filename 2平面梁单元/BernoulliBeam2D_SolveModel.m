function BernoulliBeam2D_SolveModel
%      该函数求解有限元模型，过程如下
%        1. 计算单元刚度矩阵，集成整体刚度矩阵
%        2. 计算单元的等效节点力，集成整体节点力向量
%        3. 处理约束条件，修改整体刚度矩阵和节点力向量
%        4. 求解方程组，得到整体节点位移向量

    global gNode gElement gMaterial gBC1 gNF gDF gK gDelta

    % step1. 定义整体刚度矩阵和节点力向量
    [node_number,dummy] = size( gNode ) ;
%     sparse 定义稀疏矩阵,gDelta = gK \ f也是稀疏矩阵
    gK = sparse( node_number * 3, node_number * 3 ) ;
    f = sparse( node_number * 3, 1 ) ;

    % step2. 计算单元刚度矩阵，并集成到整体刚度矩阵中
    [element_number,dummy] = size( gElement ) ;
    for ie=1:1:element_number
        k = BernoulliBeam2D_Stiffness( ie, 1 ) ;
        BernoulliBeam2D_Assembly( ie, k ) ;
    end

    % step3. 把集中力直接集成到整体节点力向量中
    [nf_number, dummy] = size( gNF ) ;
    for inf=1:1:nf_number
        n = gNF( inf, 1 ) ;
        d = gNF( inf, 2 ) ;
        f( (n-1)*3 + d ) = gNF( inf, 3 ) ;
    end

    % step4. 计算分布力的等效节点力，集成到整体节点力向量中
    [df_number, dummy] = size( gDF ) ;
    for idf = 1:1:df_number
        enf = BernoulliBeam2D_EquivalentNodeForce( gDF(idf,1), gDF(idf, 2), gDF( idf, 3), gDF( idf, 4 ) ) ;
        i = gElement( gDF(idf,1), 1 ) ;
        j = gElement( gDF(idf,1), 2 ) ;
        f( (i-1)*3+1 : (i-1)*3+3 ) = f( (i-1)*3+1 : (i-1)*3+3 ) + enf( 1:3 ) ;
        f( (j-1)*3+1 : (j-1)*3+3 ) = f( (j-1)*3+1 : (j-1)*3+3 ) + enf( 4:6 ) ;
    end
  
    % step5. 处理约束条件，修改刚度矩阵和节点力向量。采用乘大数法
    [bc_number,dummy] = size( gBC1 ) ;
    for ibc=1:1:bc_number
        n = gBC1(ibc, 1 ) ;
        d = gBC1(ibc, 2 ) ;
        m = (n-1)*3 + d ;
        f(m) = gBC1(ibc, 3)* gK(m,m) * 1e15 ;
        gK(m,m) = gK(m,m) * 1e15 ;
    end

    % step 6. 求解方程组，得到节点位移向量
    gDelta = gK \ f ;
return