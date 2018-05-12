clc
clear
Node=xlsread('Data.xls','Node');
Element=xlsread('Data.xls','Element');
Boundary=xlsread('Data.xls','Boundary');
Section=xlsread('Data.xls','Section');
Force=xlsread('Data.xls','Force');
%读取相关数据
Allstep=1000; %迭代步数
N_Node=size(Node,1); %节点个数
N_Element=size(Element,1); %单元个数
N_Force=size(Force,1); %节点力个数
N_Boundary=size(Boundary,1); %约束节点个数
Displacement=zeros(N_Node,3); %位移向量(a0)
%初始位移及转角为 0
All_Disp=zeros(N_Node,3); %初始位移和转角为零(迭代后的节点位移)
All_F=zeros(N_Node*3,1); %初始荷载向量为零 (存放节点荷载向量)
Internal_F=zeros(N_Node*3,1); %节点内力向量
ForceMatrix=zeros(N_Node*3,1); %总荷载向量
C=zeros(N_Element,2);
L=zeros(N_Element,1);
for i=1:N_Force
    ForceMatrix(Force(i,1)*3-2:Force(i,1)*3,1)=Force(i,2:4)'; %把节点荷载向量读入一个矩 阵中 , 形成列向量 =总的自由度个数
end
del=[];
for i=1:N_Boundary
    if Boundary(i,2)==0;
        m=3*Boundary(i,1)-2;





del=[del,[m]]; %受约束节点位移为 0时所对应的指标数值 1

end

if Boundary(i,3)==0;

m=3*Boundary(i,1)-1;
del=[del,[m]]; %受约束节点位移为 0时所对应的指标数值 2

end

if Boundary(i,4)==0;

m=3*Boundary(i,1);

del=[del,[m]]; %受约束节点位移为 0时所对应的指标数值 3

end

end

%求出约束节点的标号,便于刚度、荷载矩阵归 0

Update_Node=Node+Displacement(:,1:2); %更新后的节点坐标向量(x , y 坐标) Ele_F=zeros(N_Element,6); %单元节点荷载向量

ELEDISP=zeros(N_Element,6); %单元节点位移向量
Ele_F1=zeros(N_Element,6); %单元节点荷载刚度矩阵中向量

Ele1_F=zeros(1,6);

ELEDISP1=zeros(1,6);

qq(1,1)=0;

pp(1,1)=0;

for n=0:Allstep-1

n=n+1

K_Global=zeros(N_Node*3,N_Node*3); %总刚矩阵

for i=1:N_Element

N1=Element(i,1); %i节点编号

N2=Element(i,2); %j节点编号

N_Section=Element(i,3); %截面的形状控制参数

C(i,:)=Update_Node(N2,:)-Update_Node(N1,:); %a0下坐标向量增量

L(i)=norm(C(i,:)); %变形后的长度

cs=C(i,:)/L(i); %杆件的 cos 和 sin 值

E=Section(N_Section,1);

A=Section(N_Section,2);

I=Section(N_Section,3);

K_Element=beam2d_stiffness530(E,A,I,L(i),cs,Ele_F1(i,:));

K_Global=K_Global+Assemble(K_Element,Element(i,:),N_Node); %形成总刚 K0 end%整体刚度矩阵

Delta_Force=ForceMatrix/Allstep;%初始荷载向量(Q0)

Equation=[K_Global,Delta_Force];

%方程组

Disp_Transefer=ones(N_Node*3,1); %建立调节位移矩阵的位移向量

Disp_Transefer(del,:)=0; %将约束节点的位移值定为 0, 其他的定位 1
Equation(del,:)=[];%把方程中约束节点所对应的行和列去掉

Equation(:,del)=[];

%引入约束条件修改方程组

n1=size(Equation,2); % 方程组中列数

dis1=(Equation(:,1:n1-1))\Equation(:,n1); % 刚 度 矩 阵 求 逆 后 乘 以 荷 载 向 量 , Equation(:,n1)荷载向量,得到节点的位移(da0)

%求解方程组

zz=1;

%识别约束

for i=1:N_Node*3;

if Disp_Transefer(i,1)==1;
    Disp_Transefer(i,1)=dis1(zz,1); %总的位移向量

zz=zz+1;

end

end

for i=1:N_Node

Displacement(i,:)=Disp_Transefer(i*3-2:i*3,1); % 位移增量,形成 n*3的位移向量 (da0)

end

All_Disp=All_Disp+Displacement %位移向量更新得到 a1(总的位移增量) All_F=All_F+Delta_Force; %外荷载向量更新 Q1

Internal_F1=zeros(N_Node*3,1); %节点内力向量

Update_Node1=Update_Node; %C1状态

Update_Node=Node+All_Disp(:,1:2); %C2状态坐标位置更新 (改) (迭代后的位置 for i=1:N_Element

F_Global=zeros(N_Node*3,1); %全局的荷载向量

for j=1:2

ELEDISP1(j*3-2:j*3)=Displacement(Element(i,j),:); %整体 取出当前单元节点位移向量 end

N1=Element(i,1); %i节点编号

N2=Element(i,2); %j节点编号

C1(i,:)=Update_Node1(N2,:)-Update_Node1(N1,:); %a0下坐标向量增量 L1(i)=norm(C1(i,:)); %变形后的长度

cs1=C1(i,:)/L1(i); %杆件的 cos 和 sin 值

T1=[cs1(1,1),cs1(1,2),0,0,0,0;

-cs1(1,2),cs1(1,1),0,0,0,0;

0,0,1,0,0,0;

0,0,0,cs1(1,1),cs1(1,2),0;

0,0,0,-cs1(1,2),cs1(1,1),0;

0,0,0,0,0,1];

ELEDISP(i,:)=T1* ELEDISP1(:);
X1=L1(i)+ ELEDISP(i,4)- ELEDISP(i,1);

Z1= ELEDISP(i,5)-ELEDISP(i,2);

LN=(X1^2+Z1^2)^0.5;

sin1=Z1/LN;

cos1=X1/LN;

Citaloca=atan2(sin1,cos1);

Ub=LN-L1(i);
THRA=ELEDISP(i,3)-Citaloca;

THRB=ELEDISP(i,6)-Citaloca;

ELEDISP(i,:)=[0,0,THRA,Ub,0,THRB];

K_Element=beam2d_stiffness520(E,A,I,L1(i),cs1,Ele_F1(i,:)); %L(i )为 a0对应下的 Ele_F(i,:)=K_Element*ELEDISP(i,:)'; %局部坐标系下荷载(Q0)作用下的节点力 Ele_F1(i,:)= Ele_F1(i,:)+ Ele_F(i,:);

C2(i,:)=Update_Node(N2,:)-Update_Node(N1,:); %a0下坐标向量增量 L2(i)=norm(C2(i,:)); %变形后的长度

cs2=C2(i,:)/L2(i); %杆件的 cos 和 sin 值

T2=[cs2(1,1),cs2(1,2),0,0,0,0;

-cs2(1,2),cs2(1,1),0,0,0,0;

0,0,1,0,0,0;

0,0,0,cs2(1,1),cs2(1,2),0;

0,0,0,-cs2(1,2),cs2(1,1),0;

0,0,0,0,0,1];

Ele1_F(:)=T2'*Ele_F1(i,:)';%行向量

N1=Element(i,1);

N2=Element(i,2);

F_Global(3*N1-2:3*N1,1)=Ele1_F(1:3); %i节点荷载

F_Global(3*N2-2:3*N2,1)=Ele1_F(4:6); %j节点荷载

Internal_F1=Internal_F1+F_Global; %a1对应下的 P1

end

Val=Internal_F1-All_F %求出上次迭代完成时的残余应力 Q1-P1

Correc_dis=zeros(N_Node,3); %构造新的节点位移向量,每次更新

Correc_dis1=zeros(N_Node,3); %构造新的节点位移向量,叠加位移增量

N_dis=size(dis1,1); %未受约束的节点位移数目 , 不为零的节点位移数目 dis3=zeros(N_dis,1); %构造一个向量
k=0;

%修改位移矩阵形式

while norm(Val)>1e-7 & k<500 %在一个荷载增量步下进行的对此迭代

k=k+1;

K_Global=zeros(N_Node*3,N_Node*3);

for i=1:N_Element

N1=Element(i,1);

N2=Element(i,2);

N_Section=Element(i,3);

C(i,:)=Update_Node(N2,:)-Update_Node(N1,:); %a1下坐标向量增量
L(i)=norm(C(i,:)); %变形后的长度

cs=C(i,:)/L(i);

E=Section(N_Section,1);

A=Section(N_Section,2);

I=Section(N_Section,3);

K_Element=beam2d_stiffness530(E,A,I,L(i),cs,Ele_F1(i,:));

K_Global=K_Global+Assemble(K_Element,Element(i,:),N_Node); %形成总刚 k1 end

Equation=[K_Global,Val]; %在残余应力下的位移求解

Disp_Transefer=ones(N_Node*3,1);

Disp_Transefer(del,:)=0;

Equation(del,:)=[];

Equation(:,del)=[];

n1=size(Equation,2);

Dis2=-(Equation(:,1:n1-1))\Equation(:,n1) %荷位移增量 da1

zz=1;

for i=1:N_Node*3;

if Disp_Transefer(i,1)==1;

Disp_Transefer(i,1)=Dis2(zz,1);

zz=zz+1;

end

end

for i=1:N_Node

Correc_dis(i,:)=Disp_Transefer(i*3-2:i*3,1);
end

Correc_dis1= Correc_dis1+ Correc_dis;

Internal_F1=zeros(N_Node*3,1); %节点内力向量

Update_Node1=Update_Node;

Update_Node=Node+All_Disp(:,1:2)+Correc_dis1(:,1:2);%为了求切线刚度矩阵(改) a2下 if abs(Update_Node(11,2))<=1e-4&&abs(Update_Node(11,1))<=1e-4

break

end

for i=1:N_Element

F_Global=zeros(N_Node*3,1);

for j=1:2

ELEDISP1(j*3-2:j*3)=Correc_dis(Element(i,j),:); %取出当前单元节点位移向量 end
N1=Element(i,1); %i节点编号

N2=Element(i,2); %j节点编号

C1(i,:)=Update_Node1(N2,:)-Update_Node1(N1,:); %a0下坐标向量增量 L1(i)=norm(C1(i,:)); %变形后的长度

cs1=C1(i,:)/L1(i); %杆件的 cos 和 sin 值

T1=[cs1(1,1),cs1(1,2),0,0,0,0;

-cs1(1,2),cs1(1,1),0,0,0,0;

0,0,1,0,0,0;

0,0,0,cs1(1,1),cs1(1,2),0;

0,0,0,-cs1(1,2),cs1(1,1),0;

0,0,0,0,0,1];

ELEDISP(i,:)=T1* ELEDISP1(:);

X1=L1(i)+ ELEDISP(i,4)- ELEDISP(i,1);

Z1= ELEDISP(i,5)-ELEDISP(i,2);

LN=(X1^2+Z1^2)^0.5;

sin1=Z1/LN;

cos1=X1/LN;

Citaloca=atan2(sin1,cos1);

Ub=LN-L1(i);

THRA=ELEDISP(i,3)- Citaloca;

THRB=ELEDISP(i,6)- Citaloca;

ELEDISP(i,:)=[0,0,THRA,Ub,0,THRB];

K_Element=beam2d_stiffness520(E,A,I,L1(i),cs1,Ele_F1(i,:));% L(i )为 a1对应下的 Ele_F(i,:)=K_Element*ELEDISP(i,:)';

Ele_F1(i,:)= Ele_F1(i,:)+ Ele_F(i,:);

C2(i,:)=Update_Node(N2,:)-Update_Node(N1,:); %a0下坐标向量增量 L2(i)=norm(C2(i,:)); %变形后的长度
cs2=C2(i,:)/L2(i); %杆件的 cos 和 sin 值

T2=[cs2(1,1),cs2(1,2),0,0,0,0;

-cs2(1,2),cs2(1,1),0,0,0,0;

0,0,1,0,0,0;

0,0,0,cs2(1,1),cs2(1,2),0;

0,0,0,-cs2(1,2),cs2(1,1),0;

0,0,0,0,0,1];

Ele1_F(:)=T2'*Ele_F1(i,:)';%行向量

N1=Element(i,1);

N2=Element(i,2);

F_Global(3*N1-2:3*N1,1)=Ele1_F(1:3); %i节点荷载 F_Global(3*N2-2:3*N2,1)=Ele1_F(4:6); %j节点荷载 Internal_F1=Internal_F1+F_Global; %a1对应下的 P1 end

Val=Internal_F1-All_F; %荷载增量 Q1-P2

Val(del,:)=0;

end

All_Disp=All_Disp+Correc_dis1

qq(n+1,1)=All_F(N_Node*3,1);

pp(n+1,1)=All_Disp(21,2);

end

plot(1.4525,qq,'g')

text(1.3,1.5*10^4,'x=1.4525')

hold on

plot(pp,qq,'r')

title('悬臂梁加载点荷载位移曲线 ');

xlabel('位移(m ) ');

ylabel('荷载(N ) ');

legend('点 11荷载位移曲线 ');

grid


