% 主程序二:弧长法

clc

clear

Node=xlsread('Data111.xls','Node');

Element=xlsread('Data111.xls','Element');

Boundary=xlsread('Data111.xls','Boundary');
Section=xlsread('Data111.xls','Section');

Force=xlsread('Data111.xls','Force');

%读取相关数据

N_Node=size(Node,1); %节点个数

N_Element=size(Element,1); %单元个数

N_Force=size(Force,1); %节点力个数

N_Boundary=size(Boundary,1); %约束节点个数

Displacement=zeros(N_Node,3); %位移向量(a0)

%初始位移及转角为 0
All_Disp=zeros(N_Node,3); %初始位移和转角为零(迭代后的节点位移)

All_F=zeros(N_Node*3,1); %初始荷载向量为零 (存放节点荷载向量) ForceMatrix=zeros(N_Node*3,1); %总荷载向量

C=zeros(N_Element,2);

L=zeros(N_Element,1);

for i=1:N_Force

ForceMatrix(Force(i,1)*3-2:Force(i,1)*3,1)=Force(i,2:4)'; %把节点荷载向量读入一个矩 end

del=[];

for i=1:N_Boundary;

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

end %求出约束节点的标号,便于刚度、荷载矩阵归 0

Update_Node=Node+Displacement(:,1:2); %更新后的节点坐标向量(x , y 坐标) Ele_F=zeros(N_Element,6); %单元节点荷载向量

ELEDISP=zeros(N_Element,6); %单元节点位移向量

Ele_F1=zeros(N_Element,6); %单元节点荷载刚度矩阵中向量
Ele1_F=zeros(1,6); %临时存放向量

ELEDISP1=zeros(1,6); %临时存放向量

flag=1;

n=0;

DetaL=0.1; %初始迭代荷载系数

chushidian=zeros(1,N_Node*3+1);

chushidian(del)=[]; %存放 m-1点的位移 +荷载系数向量

Deta_Disp=zeros(N_Node,3); %weiyi

DetaVH=zeros(N_Node*3,1);

DetaVH(del)=[];
qq(1,1)=0;

pp(1,1)=0;

while flag==1 && n<20000

n=n+1

DetaV=zeros(N_Node*3,1);

DetaV(del)=[];

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

K_Global=K_Global+Assemble(K_Element,Element(i,:),N_Node); %形成总刚 K0 K_Global1=K_Global;

end

%整体刚度矩阵

%%%%%%%%%%%%%%%%%弧长法第一步 %%%%%%%%%%

if n==1

Equation=[K_Global1,ForceMatrix]; %方程组

Equation(del,:)=[];%把方程中约束节点所对应的行和列去掉
Equation(:,del)=[];

%引入约束条件修改方程组

n1=size(Equation,2); % 方程组中列数

dis=(Equation(:,1:n1-1))\Equation(:,n1); % 刚度矩阵求逆后乘以荷载向量 , Equation(:,n1)荷载向量,得到节点的位移(da0)

%求解方程组

Lamuda=DetaL/(sqrt(dis'*dis)+1);

dis1=Lamuda*dis;

Delta_Force=Lamuda*ForceMatrix;%初始荷载向量(Q0)

Disp_Transefer=ones(N_Node*3,1); %建立调节位移矩阵的位移向量

Disp_Transefer(del,:)=0; %将约束节点的位移值定为 0, 其他的定位 1
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

DetaV=DetaV+dis1; %用于求弧长的位移增量(第 i 荷载步的第一次迭代) DetaVH=DetaVH+dis1;

Lamuda1=Lamuda;

All_Disp=All_Disp+Displacement; %位移向量更新得到 a1(总的位移增量) All_F=All_F+Delta_Force; %外荷载向量更新 Q1

end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% if n>=2

Equation=[K_Global1,ForceMatrix];

Equation(del,:)=[];

Equation(:,del)=[];

n2=size(Equation,2); %列向量的个数

DetaVp=(Equation(:,1:n2-1))\Equation(:,n2); %更新位移

Deta_Lamuda1=-DetaL/(sqrt( DetaVp'* DetaVp+1));

Deta_Lamuda2=DetaL/(sqrt( DetaVp'* DetaVp+1));

DetaVz1=Deta_Lamuda1*DetaVp; %代表位移增量
DetaVz2=Deta_Lamuda2*DetaVp; %代表位移增量

r1=[DetaVz1',Deta_Lamuda1]'; %相对于上一个收敛点的相对向量

r2=[DetaVz2',Deta_Lamuda2]';

if fangxiangliang*r1>0

Deta_Lamuda= Deta_Lamuda1;

else

Deta_Lamuda= Deta_Lamuda2;

end
Lamuda=Lamuda+Deta_Lamuda; %更新的荷载系数

Lamuda1=Deta_Lamuda;

Deta_V=Deta_Lamuda*DetaVp; %代表位移增量

Disp_Transefer=ones(N_Node*3,1); %建立调节位移矩阵的位移向量

Disp_Transefer(del,:)=0; %将约束节点的位移值定为 0, 其他的定位 1

zz=1;

%识别约束

for i=1:N_Node*3;

if Disp_Transefer(i,1)==1;

Disp_Transefer(i,1)=Deta_V(zz,1); %总的位移向量

zz=zz+1;

end

end

for i=1:N_Node

Displacement(i,:)=Disp_Transefer(i*3-2:i*3,1); % 位移增量,形成 n*3的位移向量 (da0)

end

DetaV=DetaV+Deta_V; %用于求弧长的位移增量(第 i 荷载步的第一次迭代)

DetaVH=DetaVH+Deta_V;

All_Disp=All_Disp+Displacement; %位移向量更新得到 a1(总的位移增量) Delta_Force=Deta_Lamuda*ForceMatrix;

All_F=All_F+Delta_Force; %外荷载向量更新 Q1

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update_Node1=Update_Node; %C1状态
Update_Node=Node+All_Disp(:,1:2); %C2状态 坐标位置更新 (改) (迭代后的位置) Internal_F1=zeros(N_Node*3,1); %节点内力向量

for i=1:N_Element

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

Ele1_F(:)=T2'*Ele_F1(i,:)'; %行向量

N1=Element(i,1);

N2=Element(i,2);

F_Global(3*N1-2:3*N1,1)=Ele1_F(1:3); %i节点荷载

F_Global(3*N2-2:3*N2,1)=Ele1_F(4:6); %j节点荷载

Internal_F1=Internal_F1+F_Global; %a1对应下的 P1

end

Val=All_F-Internal_F1; %求出上次迭代完成时的残余应力 Q1-P1
Val(del,:)=0; %受约束点的残余应力项归零 , 便于下一次迭代 Correc_dis=zeros(N_Node,3); %构造新的节点位移向量,每次更新

Correc_dis1=zeros(N_Node,3); %构造新的节点位移向量,叠加位移增量 N_dis=size(dis1,1); %未受约束的节点位移数目 , 不为零的节点位移数目 dis3=zeros(N_dis,1); %构造一个向量

k=0;

%修改位移矩阵形式

% while norm(Val)/norm(Delta_Force)>=1e-7 & k<200 %在一个荷载增量步下进行 的对此迭代

while norm(Val)>=1e-4 & k<200 %在一个荷载增量步下进行的对此迭代

k=k+1;

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

K_Global=K_Global+Assemble(K_Element,Element(i,:),N_Node); %形成总刚 K0 end

Equation=[K_Global,ForceMatrix]; %求解方程组
Equation(del,:)=[];

Equation(:,del)=[];

n1=size(Equation,2); %列向量的个数

DetaV1=(Equation(:,1:n1-1))\Equation(:,n1); %zhengshu

Equation=[K_Global,Val]; %求解方程组

Equation(del,:)=[];

Equation(:,del)=[];

n2=size(Equation,2); %列向量的个数

DetaV2=(Equation(:,1:n2-1))\Equation(:,n2); %zhengshu
if k==1

Equation=[K_Global1,ForceMatrix];

Equation(del,:)=[];%把方程中约束节点所对应的行和列去掉

Equation(:,del)=[];

n1=size(Equation,2); % 方程组中列数

DetaV3=(Equation(:,1:n1-1))\Equation(:,n1); % 刚 度 矩 阵 求 逆 后 乘 以 荷 载 向 量 , Equation(:,n1)荷载向量,得到节点的位移(da0)

Deta_Lamuda=-((DetaV3')*DetaV2)/(DetaV3'*DetaV1+1)

else

Deta_Lamuda=-((DetaV')*DetaV2)/(DetaV'*DetaV1+Lamuda1)

end

DetaV3

DetaV

Lamuda=Lamuda+Deta_Lamuda; %全过程荷载因子累积

Lamuda1=Lamuda1+Deta_Lamuda; %每一荷载步下的荷载因子累积

Deta_V=Deta_Lamuda*DetaV1+DetaV2; %代表位移增量

DetaV=DetaV+Deta_V; %下一次迭代后的总位移增量

DetaVH=DetaVH+Deta_V;

All_F=Lamuda*ForceMatrix; % 更新顶部荷载

Disp_Transefer=ones(N_Node*3,1);

Disp_Transefer(del,:)=0;

zz=1;

for i=1:N_Node*3;

if Disp_Transefer(i,1)==1;

Disp_Transefer(i,1)=Deta_V(zz,1);

zz=zz+1;

end

end

for i=1:N_Node
    Correc_dis(i,:)=Disp_Transefer(i*3-2:i*3,1);

end

Correc_dis1= Correc_dis1+ Correc_dis;

Internal_F1=zeros(N_Node*3,1); %节点内力向量

Update_Node1=Update_Node;

Update_Node=Node+All_Disp(:,1:2)+Correc_dis1(:,1:2);%为 了 求 切 线 刚 度 矩 阵 (改) a2下对应的节点位置

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

F_Global(3*N1-2:3*N1,1)=Ele1_F(1:3); %i节点荷载

F_Global(3*N2-2:3*N2,1)=Ele1_F(4:6); %j节点荷载

Internal_F1=Internal_F1+F_Global; %a1对应下的 P1

end

Val=All_F-Internal_F1; %荷载增量 Q1-P2

Val(del,:)=0;

end

U0=DetaVH;

L0=Lamuda;

chushidian1=[DetaVH',Lamuda]; %存放 m 点的位移 +荷载系数向量 fangxiangliang=chushidian1-chushidian; %m-1到 m 的方向向量 chushidian=chushidian1; %更新节点向量

All_Disp=All_Disp+Correc_dis1

qq(n+1,1)= Lamuda*5000;

pp(n+1,1)=-All_Disp(9,2);

end
xlswrite('F:\data.xls',qq',pp');

plot(pp,qq,'r')

title('拱结构加载点荷载位移曲线 ');

xlabel('位移(mm ) ');

ylabel('荷载(N ) ');

legend('加载点荷载位移曲线 ');

grid
