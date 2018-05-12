clc
clear
Node=xlsread('Data.xls','Node');
Element=xlsread('Data.xls','Element');
Boundary=xlsread('Data.xls','Boundary');
Section=xlsread('Data.xls','Section');
Force=xlsread('Data.xls','Force');
%��ȡ�������
Allstep=1000; %��������
N_Node=size(Node,1); %�ڵ����
N_Element=size(Element,1); %��Ԫ����
N_Force=size(Force,1); %�ڵ�������
N_Boundary=size(Boundary,1); %Լ���ڵ����
Displacement=zeros(N_Node,3); %λ������(a0)
%��ʼλ�Ƽ�ת��Ϊ 0
All_Disp=zeros(N_Node,3); %��ʼλ�ƺ�ת��Ϊ��(������Ľڵ�λ��)
All_F=zeros(N_Node*3,1); %��ʼ��������Ϊ�� (��Žڵ��������)
Internal_F=zeros(N_Node*3,1); %�ڵ���������
ForceMatrix=zeros(N_Node*3,1); %�ܺ�������
C=zeros(N_Element,2);
L=zeros(N_Element,1);
for i=1:N_Force
    ForceMatrix(Force(i,1)*3-2:Force(i,1)*3,1)=Force(i,2:4)'; %�ѽڵ������������һ���� ���� , �γ������� =�ܵ����ɶȸ���
end
del=[];
for i=1:N_Boundary
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

end

%���Լ���ڵ�ı��,���ڸնȡ����ؾ���� 0

Update_Node=Node+Displacement(:,1:2); %���º�Ľڵ���������(x , y ����) Ele_F=zeros(N_Element,6); %��Ԫ�ڵ��������

ELEDISP=zeros(N_Element,6); %��Ԫ�ڵ�λ������
Ele_F1=zeros(N_Element,6); %��Ԫ�ڵ���ظնȾ���������

Ele1_F=zeros(1,6);

ELEDISP1=zeros(1,6);

qq(1,1)=0;

pp(1,1)=0;

for n=0:Allstep-1

n=n+1

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

K_Global=K_Global+Assemble(K_Element,Element(i,:),N_Node); %�γ��ܸ� K0 end%����նȾ���

Delta_Force=ForceMatrix/Allstep;%��ʼ��������(Q0)

Equation=[K_Global,Delta_Force];

%������

Disp_Transefer=ones(N_Node*3,1); %��������λ�ƾ����λ������

Disp_Transefer(del,:)=0; %��Լ���ڵ��λ��ֵ��Ϊ 0, �����Ķ�λ 1
Equation(del,:)=[];%�ѷ�����Լ���ڵ�����Ӧ���к���ȥ��

Equation(:,del)=[];

%����Լ�������޸ķ�����

n1=size(Equation,2); % ������������

dis1=(Equation(:,1:n1-1))\Equation(:,n1); % �� �� �� �� �� �� �� �� �� �� �� �� �� , Equation(:,n1)��������,�õ��ڵ��λ��(da0)

%��ⷽ����

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

All_Disp=All_Disp+Displacement %λ���������µõ� a1(�ܵ�λ������) All_F=All_F+Delta_Force; %������������� Q1

Internal_F1=zeros(N_Node*3,1); %�ڵ���������

Update_Node1=Update_Node; %C1״̬

Update_Node=Node+All_Disp(:,1:2); %C2״̬����λ�ø��� (��) (�������λ�� for i=1:N_Element

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

Ele1_F(:)=T2'*Ele_F1(i,:)';%������

N1=Element(i,1);

N2=Element(i,2);

F_Global(3*N1-2:3*N1,1)=Ele1_F(1:3); %i�ڵ����

F_Global(3*N2-2:3*N2,1)=Ele1_F(4:6); %j�ڵ����

Internal_F1=Internal_F1+F_Global; %a1��Ӧ�µ� P1

end

Val=Internal_F1-All_F %����ϴε������ʱ�Ĳ���Ӧ�� Q1-P1

Correc_dis=zeros(N_Node,3); %�����µĽڵ�λ������,ÿ�θ���

Correc_dis1=zeros(N_Node,3); %�����µĽڵ�λ������,����λ������

N_dis=size(dis1,1); %δ��Լ���Ľڵ�λ����Ŀ , ��Ϊ��Ľڵ�λ����Ŀ dis3=zeros(N_dis,1); %����һ������
k=0;

%�޸�λ�ƾ�����ʽ

while norm(Val)>1e-7 & k<500 %��һ�������������½��еĶԴ˵���

k=k+1;

K_Global=zeros(N_Node*3,N_Node*3);

for i=1:N_Element

N1=Element(i,1);

N2=Element(i,2);

N_Section=Element(i,3);

C(i,:)=Update_Node(N2,:)-Update_Node(N1,:); %a1��������������
L(i)=norm(C(i,:)); %���κ�ĳ���

cs=C(i,:)/L(i);

E=Section(N_Section,1);

A=Section(N_Section,2);

I=Section(N_Section,3);

K_Element=beam2d_stiffness530(E,A,I,L(i),cs,Ele_F1(i,:));

K_Global=K_Global+Assemble(K_Element,Element(i,:),N_Node); %�γ��ܸ� k1 end

Equation=[K_Global,Val]; %�ڲ���Ӧ���µ�λ�����

Disp_Transefer=ones(N_Node*3,1);

Disp_Transefer(del,:)=0;

Equation(del,:)=[];

Equation(:,del)=[];

n1=size(Equation,2);

Dis2=-(Equation(:,1:n1-1))\Equation(:,n1) %��λ������ da1

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

Internal_F1=zeros(N_Node*3,1); %�ڵ���������

Update_Node1=Update_Node;

Update_Node=Node+All_Disp(:,1:2)+Correc_dis1(:,1:2);%Ϊ�������߸նȾ���(��) a2�� if abs(Update_Node(11,2))<=1e-4&&abs(Update_Node(11,1))<=1e-4

break

end

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

F_Global(3*N1-2:3*N1,1)=Ele1_F(1:3); %i�ڵ���� F_Global(3*N2-2:3*N2,1)=Ele1_F(4:6); %j�ڵ���� Internal_F1=Internal_F1+F_Global; %a1��Ӧ�µ� P1 end

Val=Internal_F1-All_F; %�������� Q1-P2

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

title('���������ص����λ������ ');

xlabel('λ��(m ) ');

ylabel('����(N ) ');

legend('�� 11����λ������ ');

grid


