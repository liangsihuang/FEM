function UpdateInterF()
% 更新 InterF (整体坐标)
% 更新 InterF_elem (局部坐标)
global InterF InterF_elem Element
[num,~]=size(Element);
for ie=1:num
    d_inf=ElementInternalForce(ie);
    InterF_elem(:,ie)=InterF_elem(:,ie)+d_inf;
    n1=Element(ie,1);
    n2=Element(ie,2);
    T=Beam2D_TransformMatrix(ie);
    g_d_inf=T*d_inf;
    InterF(n1*3-2:n1*3)=InterF(n1*3-2:n1*3)+g_d_inf(1:3);
    InterF(n2*3-2:n2*3)=InterF(n2*3-2:n2*3)+g_d_inf(4:6);
end

end

function d_inf=ElementInternalForce(ie)
% 输入：
% ie:节点号
% 输出：
% d_inf：单元局部坐标系下的C1-C2的内力（节点力）增量
global InterF_elem
deform=RemoveRigidMotion(ie);
inf=InterF_elem(:,ie); %初应力为C1状态下的内力
T=Beam2D_TransformMatrix(ie);
inf=T'*inf; %转换到单元局部坐标
k = NonlinearBeam2D_Stiffness(ie,2,inf);%2代表局部坐标
d_inf = k * deform;
end

function deform=RemoveRigidMotion(ie)
% 输入：
% ie：单元号
% 输出：
% deform: C2局部坐标下的弹性变形 向量(6,1)
global Element dU Node
n1=Element(ie,1);
n2=Element(ie,2);
du_elem=zeros(6,1);
du_elem(1:3)=dU(n1*3-2:n1*3);
du_elem(4:6)=dU(n2*3-2:n2*3);
T=Beam2D_TransformMatrix(ie);
du_elem=T'*du_elem; %转换到单元局部坐标，C1状态
% C1到C2坐标增量各分量
ua=du_elem(1);
wa=du_elem(2);
ra=du_elem(3);
ub=du_elem(4);
wb=du_elem(5);
rb=du_elem(6);
% 求C1状态下杆件的长度, Node必须是还未叠加dU的节点坐标
x1=Node(n1,1);
y1=Node(n1,2);
x2=Node(n2,1);
y2=Node(n2,2);
L1=sqrt((x1-x2)^2+(y1-y2)^2);
% C2状态下长度，C2相对C1状态的刚体转动
dX=L1+ub-ua;
dZ=wb-wa;
L2=sqrt(dX^2+dZ^2);
r=atan(dZ/dX);
% C2局部坐标下的弹性变形
deform=zeros(6,1);
deform(4)=L2-L1;
deform(3)=ra-r;
deform(6)=rb-r;
end