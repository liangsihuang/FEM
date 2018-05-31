function SolveModel()
% 荷载增量法
for lambda=0.01:0.01:1
    NewtonRaphson(lambda);
end
% 步长太大(比如0.1)，不收敛
end


function NewtonRaphson(lambda)
global TotalU ExterF InterF dU
iter=0;
max=50;
Q=ExterF*lambda;
while(true)
    R=Q-InterF;
    R=ModifyResidual(R);
%     monitor=norm(InterF);
%     disp(monitor);
    control=norm(R);
    if (iter>max || control<1e-5)
        break;
    end
    K=GlobalStiffnessMatrix();
    [K,R]=AddBoundaryCondition(K,R,'Plain');
    dU=K\R;
%     monitor=norm(dU);
%     disp(monitor);
    UpdateInterF_elem();
    UpdateNode();
    AssembleInterF();
    TotalU=TotalU+dU;
    iter=iter+1;
end
if (iter>max)
    disp('收敛失败');
    disp(lambda);
end
if (control<1e-5)
    disp('收敛成功');
    disp(lambda);
end
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