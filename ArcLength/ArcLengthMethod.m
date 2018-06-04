function ArcLengthMethod(increNum)
global ExterF InterF dU Node lambda_total dr_total
arclength=0.05; %固定弧长
% 第一次迭代
iter=1;
Q=ExterF;
K=GlobalStiffnessMatrix();
[K,Q]=AddBoundaryCondition(K,Q,'Plain');
dU=K\Q;
dLambda=arclength/sqrt((dU'*dU+1)); %第一次迭代，根据弧长和dU反求dLambda
% lambda正负号：先假设为正，则r1与dr_total的夹角也要为正
r1=[dLambda*dU;dLambda]; %第一个弧长半径向量 ,r1=dr1
if (increNum>1)
    if r1'*dr_total<0
        dLambda=-dLambda;
        r1=-r1;
    end
end
dU=dLambda*dU; %dU要同比扩大lambda倍
du=dU; %增量步内总的位移增量，初始化，更新超平面用
lambda=dLambda; %增量步内总的lambda增量，初始化，更新超平面用
lambda_total=lambda_total+dLambda; % 求残余内力用
r=r1; %增量步内总的弧长向量，增量步内坐标（原点在m点），求r_total用
max=100;
while(true)
    UpdateStatus();
    R=Q*lambda_total-InterF; %易错，要用lambda_total
    R=ModifyResidual(R);
    converge=norm(R);
    K=GlobalStiffnessMatrix();
    [K,R]=AddBoundaryCondition(K,R,'Plain');
    temp2=K\R;
    K=GlobalStiffnessMatrix();
    [K,Q]=AddBoundaryCondition(K,Q,'Plain');
    temp3=K\Q;
    if (iter==1)
        dLambda=-dU'*temp2/(dU'*temp3+1);
    else
        dLambda=-du'*temp2/(du'*temp3+lambda);
    end
    dU=dLambda*temp3+temp2;
    du=du+dU;
    dr=[dU;dLambda];
    r=r+dr;
    lambda=lambda+dLambda;
    lambda_total=lambda_total+dLambda;
    iter=iter+1;
    if(iter>max || converge<1e-7)
        break;
    end
end
dr_total=r; %用于判断加载方向，第二个增量步开始每次第一次迭代时与r1点乘
if (converge<1e-7)
    p=['increNum=',mat2str(increNum),'收敛成功','iter=',mat2str(iter)];
    disp(p);
end
if (iter>max)
    disp('收敛失败，请手动缩小弧长');
    disp(increNum);
end
end

function UpdateStatus()
global TotalU dU
UpdateInterF_elem();
UpdateNode();
AssembleInterF();
TotalU=TotalU+dU;
end

function K=GlobalStiffnessMatrix()
global Element Node InterF_elem
[m,~] = size(Node);
K = zeros(m*3,m*3);
[num,~] = size(Element);
for ie = 1:num
    k = NonlinearBeam2D_Stiffness(ie,1,InterF_elem(:,ie));
    K = Beam2D_AssembleStiffness(K,ie,k);
end
end

function [K,R]=AddBoundaryCondition(K,R,string)
global BC1
if (strcmp(string,'Penalty'))
    [num,~] = size( BC1 ) ;
    for ibc=1:1:num
        n = BC1(ibc, 1 ) ;
        d = BC1(ibc, 2 ) ;
        m = (n-1)*3 + d ;
        R(m) = BC1(ibc, 3)* K(m,m) * 1e15 ;
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
        R = R-BC1(ibc,3)*K(:,m);
        R(m) = BC1(ibc, 3);
        K(:,m) = zeros(p,1);
        K(m,:) = zeros(1,p);
        K(m,m) = 1.0;
    end
end
end


function R=ModifyResidual(R)
global BC1
[num,~] = size( BC1 ) ;
for i=1:1:num
    n = BC1(i, 1);
    d = BC1(i, 2);
    m = (n-1)*3 + d;
    R(m) = BC1(i, 3);
end
end

