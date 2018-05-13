function NewtonRaphson(Q)
conv_tol=1*10^(-5);
iter_max=20;
conv=1000;
iter=1;
R=Q;
Kt=TangentStiffness();
while conv>conv_tol && iter<iter_max
    del_u=Kt\R;
end

end

function Kt=TangentStiffness()
global gNode gElement
[node_number,~]=size(gNode);
Kt=sparse(node_number*3,node_number*3);
[element_number,~]=size(gElement);
for ie=1:1:element_number
    k = NonlinearBeam2D_Stiffness( ie, 1 ) ;
    Kt=Beam2D_Assembly(ie,k,Kt);
end
end

function F=UpdatedInternalForce(del_u,F)
global gElement
[element_number,~]=size(gElement);
for ie=1:1:element_number
    i = gElement( ie, 1 ) ;
    j = gElement( ie, 2 ) ;
% 提取单元位移增量，注意：位移增量是整体坐标！
    del_ue = zeros( 6, 1 );
    del_ue( 1:3 ) = del_u( (i-1)*3+1:(i-1)*3+3 );
    del_ue( 4:6 ) = del_u( (j-1)*3+1:(j-1)*3+3 );
% 去除刚体位移
    del_ue=RemoveRigidMotion(ie,del_ue);
    k=NonlinearBeam2D_Stiffness( ie, 1 );
    del_F=k*del_ue;
% 组装到整体内力矩阵
    F((i-1)*3+1:(i-1)*3+3)=F((i-1)*3+1:(i-1)*3+3)+del_F(1:3);
    F((j-1)*3+1:(j-1)*3+3)=F((j-1)*3+1:(j-1)*3+3)+del_F(4:6);
end
end

function del_ue=RemoveRigidMotion(ie,del_ue)
global gNode
xi = gNode( gElement( ie, 1 ), 1 ) ;
yi = gNode( gElement( ie, 1 ), 2 ) ;
xj = gNode( gElement( ie, 2 ), 1 ) ;
yj = gNode( gElement( ie, 2 ), 2 ) ;
L_old = ( (xj-xi)^2 + (yj-yi)^2 )^(1/2) ;
ui=del_ue(1);wi=del_ue(2);thetai=del_ue(3);
uj=del_ue(4);wj=del_ue(5);thetaj=del_ue(6);
X=L_old+uj-ui;
Z=wj-wi;
L_new=sqrt(X^2+Z^2);
elongation=L_new-L_old;
rotation_rigid=atan(Z/X);
rotation_i=thetai-rotation_rigid;
rotation_j=thetaj-rotation_rigid;
del_ue(1)=0;del_ue(2)=0;del_ue(3)=rotation_i;
del_ue(4)=elongation;del_ue(5)=0;del_ue(6)=rotation_j;
end

