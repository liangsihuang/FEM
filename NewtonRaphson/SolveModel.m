function SolveModel()
for lambda=0:0.1:1
    NewtonRaphson(lambda);
end
end

function AddBoundaryCondition()
end


function NewtonRaphson(lambda)
global Disp Node
Q=NF*lambda;

R=Q-F;
K=GlobalStiffnessMatrix();
AddBoundaryCondition();
du=K\R;
Node=Node+du(:,1:2);
Disp=Disp+du;
F=InternalForce();
end

function K=GlobalStiffnessMatrix()
[element_number,~] = size( gElement ) ;
    for ie=1:1:element_number
        k = BernoulliBeam2D_Stiffness( ie, 1 ) ;
        BernoulliBeam2D_Assembly( ie, k ) ;
    end
end