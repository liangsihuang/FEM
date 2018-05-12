% �������:������

clc

clear

Node=xlsread('Data111.xls','Node');

Element=xlsread('Data111.xls','Element');

Boundary=xlsread('Data111.xls','Boundary');
Section=xlsread('Data111.xls','Section');

Force=xlsread('Data111.xls','Force');

%��ȡ�������

N_Node=size(Node,1); %�ڵ����

N_Element=size(Element,1); %��Ԫ����

N_Force=size(Force,1); %�ڵ�������

N_Boundary=size(Boundary,1); %Լ���ڵ����

Displacement=zeros(N_Node,3); %λ������(a0)

%��ʼλ�Ƽ�ת��Ϊ 0
All_Disp=zeros(N_Node,3); %��ʼλ�ƺ�ת��Ϊ��(������Ľڵ�λ��)

All_F=zeros(N_Node*3,1); %��ʼ��������Ϊ�� (��Žڵ��������) ForceMatrix=zeros(N_Node*3,1); %�ܺ�������

C=zeros(N_Element,2);

L=zeros(N_Element,1);

for i=1:N_Force

ForceMatrix(Force(i,1)*3-2:Force(i,1)*3,1)=Force(i,2:4)'; %�ѽڵ������������һ���� end

del=[];

for i=1:N_Boundary;

if Boundary(i,2)==0;

m=3*Boundary(i,1)-2;

del=[del,[m]]; %��Լ���ڵ�λ��Ϊ 0ʱ����Ӧ��ָ����ֵ 1

end

if Boundary(i,3)==0;

m=3*Boundary(i,1)-1;

del=[del,[m]]; %��Լ���ڵ�λ��Ϊ 0ʱ����Ӧ��ָ����ֵ 2

end

if Boundary(i,4)==0;

m=3*Boundary(i,1);

del=[del,[m]]; %��Լ���ڵ�λ��Ϊ 0ʱ����Ӧ��ָ����ֵ 3

end

end %���Լ���ڵ�ı��,���ڸնȡ����ؾ���� 0

Update_Node=Node+Displacement(:,1:2); %���º�Ľڵ���������(x , y ����) Ele_F=zeros(N_Element,6); %��Ԫ�ڵ��������

ELEDISP=zeros(N_Element,6); %��Ԫ�ڵ�λ������

Ele_F1=zeros(N_Element,6); %��Ԫ�ڵ���ظնȾ���������
Ele1_F=zeros(1,6); %��ʱ�������

ELEDISP1=zeros(1,6); %��ʱ�������

flag=1;

n=0;

DetaL=0.1; %��ʼ��������ϵ��

chushidian=zeros(1,N_Node*3+1);

chushidian(del)=[]; %��� m-1���λ�� +����ϵ������

Deta_Disp=zeros(N_Node,3); %weiyi

DetaVH=zeros(N_Node*3,1);

DetaVH(del)=[];
qq(1,1)=0;

pp(1,1)=0;

while flag==1 && n<20000

n=n+1

DetaV=zeros(N_Node*3,1);

DetaV(del)=[];

K_Global=zeros(N_Node*3,N_Node*3); %�ܸվ���

for i=1:N_Element

N1=Element(i,1); %i�ڵ���

N2=Element(i,2); %j�ڵ���

N_Section=Element(i,3); %�������״���Ʋ���

C(i,:)=Update_Node(N2,:)-Update_Node(N1,:); %a0��������������

L(i)=norm(C(i,:)); %���κ�ĳ���

cs=C(i,:)/L(i); %�˼��� cos �� sin ֵ

E=Section(N_Section,1);

A=Section(N_Section,2);

I=Section(N_Section,3);

K_Element=beam2d_stiffness530(E,A,I,L(i),cs,Ele_F1(i,:));

K_Global=K_Global+Assemble(K_Element,Element(i,:),N_Node); %�γ��ܸ� K0 K_Global1=K_Global;

end

%����նȾ���

%%%%%%%%%%%%%%%%%��������һ�� %%%%%%%%%%

if n==1

Equation=[K_Global1,ForceMatrix]; %������

Equation(del,:)=[];%�ѷ�����Լ���ڵ�����Ӧ���к���ȥ��
Equation(:,del)=[];

%����Լ�������޸ķ�����

n1=size(Equation,2); % ������������

dis=(Equation(:,1:n1-1))\Equation(:,n1); % �նȾ����������Ժ������� , Equation(:,n1)��������,�õ��ڵ��λ��(da0)

%��ⷽ����

Lamuda=DetaL/(sqrt(dis'*dis)+1);

dis1=Lamuda*dis;

Delta_Force=Lamuda*ForceMatrix;%��ʼ��������(Q0)

Disp_Transefer=ones(N_Node*3,1); %��������λ�ƾ����λ������

Disp_Transefer(del,:)=0; %��Լ���ڵ��λ��ֵ��Ϊ 0, �����Ķ�λ 1
zz=1;

%ʶ��Լ��

for i=1:N_Node*3;

if Disp_Transefer(i,1)==1;

Disp_Transefer(i,1)=dis1(zz,1); %�ܵ�λ������

zz=zz+1;

end

end

for i=1:N_Node

Displacement(i,:)=Disp_Transefer(i*3-2:i*3,1); % λ������,�γ� n*3��λ������ (da0)

end

DetaV=DetaV+dis1; %�����󻡳���λ������(�� i ���ز��ĵ�һ�ε���) DetaVH=DetaVH+dis1;

Lamuda1=Lamuda;

All_Disp=All_Disp+Displacement; %λ���������µõ� a1(�ܵ�λ������) All_F=All_F+Delta_Force; %������������� Q1

end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% if n>=2

Equation=[K_Global1,ForceMatrix];

Equation(del,:)=[];

Equation(:,del)=[];

n2=size(Equation,2); %�������ĸ���

DetaVp=(Equation(:,1:n2-1))\Equation(:,n2); %����λ��

Deta_Lamuda1=-DetaL/(sqrt( DetaVp'* DetaVp+1));

Deta_Lamuda2=DetaL/(sqrt( DetaVp'* DetaVp+1));

DetaVz1=Deta_Lamuda1*DetaVp; %����λ������
DetaVz2=Deta_Lamuda2*DetaVp; %����λ������

r1=[DetaVz1',Deta_Lamuda1]'; %�������һ����������������

r2=[DetaVz2',Deta_Lamuda2]';

if fangxiangliang*r1>0

Deta_Lamuda= Deta_Lamuda1;

else

Deta_Lamuda= Deta_Lamuda2;

end
Lamuda=Lamuda+Deta_Lamuda; %���µĺ���ϵ��

Lamuda1=Deta_Lamuda;

Deta_V=Deta_Lamuda*DetaVp; %����λ������

Disp_Transefer=ones(N_Node*3,1); %��������λ�ƾ����λ������

Disp_Transefer(del,:)=0; %��Լ���ڵ��λ��ֵ��Ϊ 0, �����Ķ�λ 1

zz=1;

%ʶ��Լ��

for i=1:N_Node*3;

if Disp_Transefer(i,1)==1;

Disp_Transefer(i,1)=Deta_V(zz,1); %�ܵ�λ������

zz=zz+1;

end

end

for i=1:N_Node

Displacement(i,:)=Disp_Transefer(i*3-2:i*3,1); % λ������,�γ� n*3��λ������ (da0)

end

DetaV=DetaV+Deta_V; %�����󻡳���λ������(�� i ���ز��ĵ�һ�ε���)

DetaVH=DetaVH+Deta_V;

All_Disp=All_Disp+Displacement; %λ���������µõ� a1(�ܵ�λ������) Delta_Force=Deta_Lamuda*ForceMatrix;

All_F=All_F+Delta_Force; %������������� Q1

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update_Node1=Update_Node; %C1״̬
Update_Node=Node+All_Disp(:,1:2); %C2״̬ ����λ�ø��� (��) (�������λ��) Internal_F1=zeros(N_Node*3,1); %�ڵ���������

for i=1:N_Element

F_Global=zeros(N_Node*3,1); %ȫ�ֵĺ�������

for j=1:2

ELEDISP1(j*3-2:j*3)=Displacement(Element(i,j),:); %���� ȡ����ǰ��Ԫ�ڵ�λ������ end

N1=Element(i,1); %i�ڵ���

N2=Element(i,2); %j�ڵ���

C1(i,:)=Update_Node1(N2,:)-Update_Node1(N1,:); %a0�������������� L1(i)=norm(C1(i,:)); %���κ�ĳ���

cs1=C1(i,:)/L1(i); %�˼��� cos �� sin ֵ
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

K_Element=beam2d_stiffness520(E,A,I,L1(i),cs1,Ele_F1(i,:)); %L(i )Ϊ a0��Ӧ�µ� Ele_F(i,:)=K_Element*ELEDISP(i,:)'; %�ֲ�����ϵ�º���(Q0)�����µĽڵ��� Ele_F1(i,:)= Ele_F1(i,:)+ Ele_F(i,:);

C2(i,:)=Update_Node(N2,:)-Update_Node(N1,:); %a0�������������� L2(i)=norm(C2(i,:)); %���κ�ĳ���

cs2=C2(i,:)/L2(i); %�˼��� cos �� sin ֵ

T2=[cs2(1,1),cs2(1,2),0,0,0,0;

-cs2(1,2),cs2(1,1),0,0,0,0;

0,0,1,0,0,0;

0,0,0,cs2(1,1),cs2(1,2),0;
0,0,0,-cs2(1,2),cs2(1,1),0;

0,0,0,0,0,1];

Ele1_F(:)=T2'*Ele_F1(i,:)'; %������

N1=Element(i,1);

N2=Element(i,2);

F_Global(3*N1-2:3*N1,1)=Ele1_F(1:3); %i�ڵ����

F_Global(3*N2-2:3*N2,1)=Ele1_F(4:6); %j�ڵ����

Internal_F1=Internal_F1+F_Global; %a1��Ӧ�µ� P1

end

Val=All_F-Internal_F1; %����ϴε������ʱ�Ĳ���Ӧ�� Q1-P1
Val(del,:)=0; %��Լ����Ĳ���Ӧ������� , ������һ�ε��� Correc_dis=zeros(N_Node,3); %�����µĽڵ�λ������,ÿ�θ���

Correc_dis1=zeros(N_Node,3); %�����µĽڵ�λ������,����λ������ N_dis=size(dis1,1); %δ��Լ���Ľڵ�λ����Ŀ , ��Ϊ��Ľڵ�λ����Ŀ dis3=zeros(N_dis,1); %����һ������

k=0;

%�޸�λ�ƾ�����ʽ

% while norm(Val)/norm(Delta_Force)>=1e-7 & k<200 %��һ�������������½��� �ĶԴ˵���

while norm(Val)>=1e-4 & k<200 %��һ�������������½��еĶԴ˵���

k=k+1;

K_Global=zeros(N_Node*3,N_Node*3); %�ܸվ���

for i=1:N_Element

N1=Element(i,1); %i�ڵ���

N2=Element(i,2); %j�ڵ���

N_Section=Element(i,3); %�������״���Ʋ���

C(i,:)=Update_Node(N2,:)-Update_Node(N1,:); %a0��������������

L(i)=norm(C(i,:)); %���κ�ĳ���

cs=C(i,:)/L(i); %�˼��� cos �� sin ֵ

E=Section(N_Section,1);

A=Section(N_Section,2);

I=Section(N_Section,3);

K_Element=beam2d_stiffness530(E,A,I,L(i),cs,Ele_F1(i,:));

K_Global=K_Global+Assemble(K_Element,Element(i,:),N_Node); %�γ��ܸ� K0 end

Equation=[K_Global,ForceMatrix]; %��ⷽ����
Equation(del,:)=[];

Equation(:,del)=[];

n1=size(Equation,2); %�������ĸ���

DetaV1=(Equation(:,1:n1-1))\Equation(:,n1); %zhengshu

Equation=[K_Global,Val]; %��ⷽ����

Equation(del,:)=[];

Equation(:,del)=[];

n2=size(Equation,2); %�������ĸ���

DetaV2=(Equation(:,1:n2-1))\Equation(:,n2); %zhengshu
if k==1

Equation=[K_Global1,ForceMatrix];

Equation(del,:)=[];%�ѷ�����Լ���ڵ�����Ӧ���к���ȥ��

Equation(:,del)=[];

n1=size(Equation,2); % ������������

DetaV3=(Equation(:,1:n1-1))\Equation(:,n1); % �� �� �� �� �� �� �� �� �� �� �� �� �� , Equation(:,n1)��������,�õ��ڵ��λ��(da0)

Deta_Lamuda=-((DetaV3')*DetaV2)/(DetaV3'*DetaV1+1)

else

Deta_Lamuda=-((DetaV')*DetaV2)/(DetaV'*DetaV1+Lamuda1)

end

DetaV3

DetaV

Lamuda=Lamuda+Deta_Lamuda; %ȫ���̺��������ۻ�

Lamuda1=Lamuda1+Deta_Lamuda; %ÿһ���ز��µĺ��������ۻ�

Deta_V=Deta_Lamuda*DetaV1+DetaV2; %����λ������

DetaV=DetaV+Deta_V; %��һ�ε��������λ������

DetaVH=DetaVH+Deta_V;

All_F=Lamuda*ForceMatrix; % ���¶�������

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

Internal_F1=zeros(N_Node*3,1); %�ڵ���������

Update_Node1=Update_Node;

Update_Node=Node+All_Disp(:,1:2)+Correc_dis1(:,1:2);%Ϊ �� �� �� �� �� �� �� �� (��) a2�¶�Ӧ�Ľڵ�λ��

for i=1:N_Element
    F_Global=zeros(N_Node*3,1);

for j=1:2

ELEDISP1(j*3-2:j*3)=Correc_dis(Element(i,j),:); %ȡ����ǰ��Ԫ�ڵ�λ������ end

N1=Element(i,1); %i�ڵ���

N2=Element(i,2); %j�ڵ���

C1(i,:)=Update_Node1(N2,:)-Update_Node1(N1,:); %a0�������������� L1(i)=norm(C1(i,:)); %���κ�ĳ���

cs1=C1(i,:)/L1(i); %�˼��� cos �� sin ֵ

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

K_Element=beam2d_stiffness520(E,A,I,L1(i),cs1,Ele_F1(i,:));% L(i )Ϊ a1��Ӧ�µ� Ele_F(i,:)=K_Element*ELEDISP(i,:)';

Ele_F1(i,:)= Ele_F1(i,:)+ Ele_F(i,:);

C2(i,:)=Update_Node(N2,:)-Update_Node(N1,:); %a0�������������� L2(i)=norm(C2(i,:)); %���κ�ĳ���

cs2=C2(i,:)/L2(i); %�˼��� cos �� sin ֵ

T2=[cs2(1,1),cs2(1,2),0,0,0,0;

-cs2(1,2),cs2(1,1),0,0,0,0;

0,0,1,0,0,0;

0,0,0,cs2(1,1),cs2(1,2),0;

0,0,0,-cs2(1,2),cs2(1,1),0;
0,0,0,0,0,1];

Ele1_F(:)=T2'*Ele_F1(i,:)';%������

N1=Element(i,1);

N2=Element(i,2);

F_Global(3*N1-2:3*N1,1)=Ele1_F(1:3); %i�ڵ����

F_Global(3*N2-2:3*N2,1)=Ele1_F(4:6); %j�ڵ����

Internal_F1=Internal_F1+F_Global; %a1��Ӧ�µ� P1

end

Val=All_F-Internal_F1; %�������� Q1-P2

Val(del,:)=0;

end

U0=DetaVH;

L0=Lamuda;

chushidian1=[DetaVH',Lamuda]; %��� m ���λ�� +����ϵ������ fangxiangliang=chushidian1-chushidian; %m-1�� m �ķ������� chushidian=chushidian1; %���½ڵ�����

All_Disp=All_Disp+Correc_dis1

qq(n+1,1)= Lamuda*5000;

pp(n+1,1)=-All_Disp(9,2);

end
xlswrite('F:\data.xls',qq',pp');

plot(pp,qq,'r')

title('���ṹ���ص����λ������ ');

xlabel('λ��(mm ) ');

ylabel('����(N ) ');

legend('���ص����λ������ ');

grid
