function enf = NonlinearBeam2D_ElementNodeForce( ie )
% 输入:
% ie:节点号
% 输出：
% enf：单元局部坐标系下的节点力
deform=RemoveRigidMotion(ie);
inf=zeros(6,1); %假设初应力为0
k = NonlinearBeam2D_Stiffness(ie,2,inf);%2代表局部坐标
enf = k * deform;
return