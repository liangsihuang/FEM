function [lambda, r_total]=ArcLengthMethod(increNum, lambda, r_total)
arclength=0.05; %��ʼ������
global ExterF InterF dU
iter=0;
max=50;
Q=ExterF;
% ��һ�ε���
K=GlobalStiffnessMatrix();
[K,Q]=AddBoundaryCondition(K,Q,'Plain');
dU=K\Q;
temp1=dU;
lambda=arclength/sqrt((dU'*dU+1)); %��һ�ε��������ݻ�����dU����lambda
dU=lambda*dU; %dUҪͬ������lambda��
r1=[dU,lambda]; %��һ�������뾶���� ,r1=dr1
if (increNum>1)
    if r1'*r_total<0
        lambda=-lambda;
        dU=-dU;
        r1=-r1;
        Q=-Q;
        temp1=-temp1;
    end
end
r_total=r_total+r1;
UpdateStatus();
R=Q*lambda-InterF;
R=ModifyResidual(R);
converge=norm(R);
while(iter>max || converge<1e-7)
    
    K=GlobalStiffnessMatrix();
    [K,R]=AddBoundaryCondition(K,R,'Plain');
    temp2=K\R;
    K=GlobalStiffnessMatrix();
    [K,Q]=AddBoundaryCondition(K,Q,'Plain');
    temp3=K\Q;
    dLambda=-temp1'*temp2/(temp1'*temp3+1); 
    dU=dLambda*temp3+temp2;
    dr=[dU,dLambda];
    r_total=r_total+dr;
    lambda=lambda+dLambda;
    iter=iter+1;
end
if (converge<1e-7)
    disp('�����ɹ�');
    disp(lambda);
end
if (iter>max)
    disp('����ʧ�ܣ����ֶ���С����');
    disp(arclength);
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
if (strcmp(string,'Plain')) %���л��з���0-1��
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

