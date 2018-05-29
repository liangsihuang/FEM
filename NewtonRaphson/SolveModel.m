function SolveModel()
% for lambda=0.1:0.1:1
lambda=0.1;
    NewtonRaphson(lambda);
% end
end


function NewtonRaphson(lambda)
global Disp ExterF InterF dU
iter=0;
max=50;
Q=ExterF*lambda;
while(true)
    R=Q-InterF;
    control=norm(R);
    if (iter>max || control<1e-5)
        break;
    end
    K=GlobalStiffnessMatrix();
    [K,R]=AddBoundaryCondition(K,R);
    dU=K\R;
    UpdateInterF();
    UpdateNode();
    Disp=Disp+dU;
    iter=iter+1;
end
if (iter>max)
    disp(' ’¡≤ ß∞‹');
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

function [K,R]=AddBoundaryCondition(K,R)
global BC1
[num,~] = size( BC1 ) ;
for i=1:1:num
    n = BC1(i, 1 );
    d = BC1(i, 2 );
    m = (n-1)*3 + d ;
    R(m) = BC1(i, 3)* K(m,m) * 1e15 ;
    K(m,m) = K(m,m) * 1e15 ;
end
end