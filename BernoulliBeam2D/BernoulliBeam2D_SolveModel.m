function BernoulliBeam2D_SolveModel
%      该函数求解有限元模型，过程如下
%        1. 计算单元刚度矩阵，集成整体刚度矩阵
%        2. 计算单元的等效节点力，集成整体节点力向量
%        3. 处理约束条件，修改整体刚度矩阵和节点力向量
%        4. 求解方程组，得到整体节点位移向量

global Node Element BC1 NF DF K dU f_onlyNF

% step1. 定义整体刚度矩阵和节点力向量
[node_number,~] = size( Node ) ;
%     sparse 定义稀疏矩阵,gDelta = gK \ f也是稀疏矩阵
%     K = sparse( node_number * 3, node_number * 3 ) ;
%     f = sparse( node_number * 3, 1 ) ;
K = zeros( node_number * 3, node_number * 3 ) ;
f = zeros( node_number * 3, 1 ) ;
    
% step2. 计算单元刚度矩阵，并集成到整体刚度矩阵中
GlobalStiffnessMatrix()

% step3. 把集中力直接集成到整体节点力向量中
[nf_number, ~] = size( NF ) ;
for inf=1:1:nf_number
    n = NF( inf, 1 ) ;
    d = NF( inf, 2 ) ;
    f( (n-1)*3 + d ) = NF( inf, 3 ) ;
end

f_onlyNF=f;
    
% step4. 计算分布力的等效节点力，集成到整体节点力向量中
[df_number, ~] = size( DF ) ;
for idf = 1:1:df_number
    enf = BernoulliBeam2D_EquivalentNodeForce( DF(idf,1), DF(idf, 2), DF( idf, 3), DF( idf, 4));
    i = Element( DF(idf,1), 1 ) ;
    j = Element( DF(idf,1), 2 ) ;
    f( (i-1)*3+1 : (i-1)*3+3 ) = f( (i-1)*3+1 : (i-1)*3+3 ) + enf( 1:3 ) ;
    f( (j-1)*3+1 : (j-1)*3+3 ) = f( (j-1)*3+1 : (j-1)*3+3 ) + enf( 4:6 ) ;
end
  
% step5. 处理约束条件，修改刚度矩阵和节点力向量。采用乘大数法
[K,f]=AddBoundaryCondition(K,f,'Plain');

% step 6. 求解方程组，得到节点位移向量
dU = K \ f ;
end

function GlobalStiffnessMatrix()
global Element
[num,~] = size( Element ) ;
for ie=1:num
    ibc = BernoulliBeam2D_Stiffness( ie, 1 ) ;
    BernoulliBeam2D_Assembly( ie, ibc ) ;
end
end

function [K,f]=AddBoundaryCondition(K,f,string)
global BC1
if (strcmp(string,'Penalty'))
    [num,~] = size( BC1 ) ;
    for ibc=1:1:num
        n = BC1(ibc, 1 ) ;
        d = BC1(ibc, 2 ) ;
        m = (n-1)*3 + d ;
        f(m) = BC1(ibc, 3)* K(m,m) * 1e15 ;
        K(m,m) = K(m,m) * 1e15 ;
    end
end
if (strcmp(string,'Plain')) %划行划列法，0-1法
    [num,~] = size( BC1 ) ;
    [p,~]=size(K);
    for ibc=1:1:num
        n = BC1(ibc, 1 );
        d = BC1(ibc, 2 );
        m = (n-1)*3 + d;
        f = f-BC1(ibc,3)*K(:,m);
        f(m) = BC1(ibc, 3);
        K(:,m) = zeros(p,1);
        K(m,:) = zeros(1,p);
        K(m,m) = 1.0;
    end
end
end