% Newton_Raphson�� ������ 
clear;
clc;
% ��������
Node=xlsread('Data.xlsx','Node'); % ����ڵ�������Ϣ
Element=xlsread('Data.xlsx','Element'); % ���뵥Ԫ��ż���Ԫ����ڵ���
Boundary=xlsread('Data.xlsx','Boundary'); % ����߽������������н�1�ڵ�̶�Լ����
Property=xlsread('Data.xlsx','Property'); % �����������һ��Ϊ����ģ�����ڶ���Ϊ���������������Ϊ������Ծ�
External_Load=xlsread('Data.xlsx','External_Load'); % ��������أ��γɺ�����Ϣ���󣨱�������������������˲�������أ�

Num_node=size(Node,1); % �ڵ�����
Num_element=size(Element,1); % ��Ԫ����
Num_el=size(External_Load,1); % �ⲿ���ú�������
Num_bn=size(Boundary,1); % ��Լ���ڵ�����
Displacement=zeros(Num_node,3); % ����λ�ƴ洢���󣬳�ʼλ�ƾ�Ϊ0
Internal_force=zeros(Num_node*3,1); % ���ɿյĽڵ������洢����
Dis=zeros(Num_node,3); % ���ɿյĽڵ�λ��ת�ǵ����洢����
Ex_Load=zeros(Num_node*3,1); % ���ɿյĳ�ʼ�ڵ�����ش洢����
Total_load=zeros(Num_node*3,1); % �ܺ�������
Total_step=1000; % �����ܴ���

for i=1:Num_el
    Total_load(External_Load(i,1)*3-2:External_Load(i,1)*3,1)=External_Load(i,2:4)'; % ��������Ϣ����ת��Ϊ���е�����ؾ���
end

% �߽�����������Ϣ�����Լ���ڵ�ı�ţ����ڸնȡ����ؾ����0
correct=[]; % ���ɿվ���
for i=1:Num_bn
    if Boundary(i,2)==0
       m=Boundary(i,1)*3-2;    % ��Լ���ڵ�x��λ��Ϊ0ʱ����Ӧ�����б��Ϊ 3*�ڵ���-2
       correct=[correct,[m]];  % ��¼���б��
    end
    if Boundary(i,3)==0
       m=Boundary(i,1)*3-1;    % ��Լ���ڵ�y��λ��Ϊ0ʱ����Ӧ�����б��Ϊ 3*�ڵ���-1
       correct=[correct,[m]];  % ��¼���б�� 
    end
    if Boundary(i,4)==0
       m=Boundary(i,1)*3;      % ��Լ���ڵ�ת��Ϊ0ʱ����Ӧ�����б��Ϊ 3*�ڵ���
       correct=[correct,[m]];  % ��¼���б��
    end
end

Update_Node=Node+Dis(:,1:2); % ����λ����Ϣ���µĽڵ�����,x������y����
Coordinate=zeros(Num_node,2); % ��Ԫ�ڵ����������洢����
Length1=zeros(Num_node,1); % ��Ԫ���ȴ���
Element_Load=zeros(Num_element,6); % ��Ԫ�ڵ������
Element_Disp=zeros(Num_element,6); % ��Ԫ�ڵ�λ��
Element_Load_Geo=zeros(Num_element,6); % ��Ԫ�ڵ���������뼸�ξ����йص������
Element_Load_2=zeros(1,6); % ������Ԫ�ڵ�����ʱ�洢����
Ele_dis=zeros(1,6); % ������Ԫ�ڵ�λ����ʱ�洢����
x(1,1)=0; % ��ͼ����
y(1,1)=0; % ��ͼ����

for n=1:Total_step
    n=n+1;
    K=zeros(Num_node*3,Num_node*3); % �γɿյ��ܸվ���
    for i=1:Num_element    % �����ܸ� 
        N1=Element(i,1);   % ��ȡ��Ԫi�ڵ���
        N2=Element(i,2);   % ��ȡ��Ԫj�ڵ���
        N_E=Element(i,3);  % ��ȡ��Ԫ���
        E=Property(N_E,1); % ��ȡ��Ԫ���ϵ���ģ������λΪN/m^2
        A=Property(N_E,2); % ��ȡ��Ԫ�����������λΪm^2
        I=Property(N_E,3); % ��ȡ��Ԫ���Ծأ���λΪm^4
        Coordinate(i,:)=Update_Node(N2,:)-Update_Node(N1,:); % �������ɵĵ�Ԫ��������
        Length1(i)=norm(Coordinate(i,:)); % �������ɵĵ�Ԫ���ȱ仯
        cs=Coordinate(i,:)/Length1(i); % ��Ԫת�ǵ�cos��sinֵ
        KE=beam2d_stiffness(E,A,I,Length1(i),cs,Element_Load_Geo(i,:));
        K=K+Assemble(KE,Element(i,:),Num_node); % �õ��ܸ�
    end
    Load_step=Total_load/Total_step;  % ������س�ʼ��������С
    Equation_save=[K,Load_step];      % �洢��������ϵ��
    Disp_correct=ones(Num_node*3,1);  % ����λ�ƾ������õ�λ������
    Disp_correct(correct,:)=0;        % ����Լ����λ�Ƶ�ֵΪ0�������Ϊ1
    % ����߽�������������
    Equation_save(correct,:)=[];      % ���������ܵ�Լ����λ�ƶ�Ӧ���л�ȥ
    Equation_save(:,correct)=[];      % ���������ܵ�Լ����λ�ƶ�Ӧ���л�ȥ
    Num_column=size(Equation_save,2); % ��ȡ����ϵ���е�����
    % ��ⷽ���飬�õ��ڵ�λ��
    dis=Equation_save(:,1:Num_column-1)\Equation_save(:,Num_column);
    aa=1;
    for i=1:Num_node*3 % ��λ����ȡ��д��λ���������У��γ��������λ��������
        if Disp_correct(i,1)==1
        	 Disp_correct(i,1)=dis(aa,1);
        	 aa=aa+1;
        end	 
    end
    for i=1:Num_node
        Dis(i,:)=Disp_correct(3*i-2:3*i,1); % ���������λ��д��λ�ƾ���
    end
    Displacement=Dis+Displacement; % ����λ�ƾ��󣬵õ���λ������
    Ex_Load=Ex_Load+Load_step; % �������������
    Internal_force=zeros(Num_node*3,1); % ���ɿյĽڵ������洢����
    Update_Node_1=Update_Node; % C1״̬���ò�����ǰ��λ��
    Update_Node=Node+Displacement(:,1:2); % C2״̬���ò��������λ��
    for i=1:Num_element
        Ex_load_global=zeros(Num_node*3,1); % ����ȫ�ֽڵ��������
        for j=1:2
            % ������λ�ƾ����еõ���ǰ��Ԫ�Ľڵ�λ������
            Ele_dis(j*3-2:j*3)=Dis(Element(i,j),:);
        end
        N1=Element(i,1);   % ��ȡ��Ԫi�ڵ���
        N2=Element(i,2);   % ��ȡ��Ԫj�ڵ���
        % C1״̬
        Coordinate_1(i,:)=Update_Node_1(N2,:)-Update_Node_1(N1,:); % ��Ԫ����ǰ����
        Length_1(i)=norm(Coordinate_1(i,:));  % ��Ԫ����ǰ�ĳ���
        cs_1=Coordinate_1(i,:)/Length_1(i); % ��Ԫת�ǵ�cos��sinֵ
        T_1=[cs_1(1,1),cs_1(1,2),0,0,0,0;...
            -cs_1(1,2),cs_1(1,1),0,0,0,0;...
            0,0,1,0,0,0;...
            0,0,0,cs_1(1,1),cs_1(1,2),0;...
            0,0,0,-cs_1(1,2),cs_1(1,1),0;...
            0,0,0,0,0,1];
        Element_Disp(i,:)=T_1*Ele_dis(:);
        X_1=Length_1(i)+Element_Disp(i,4)-Element_Disp(i,1);
        Z_1=Element_Disp(i,5)-Element_Disp(i,2);
        L_1=(X_1^2+Z_1^2)^0.5; % ��Ԫ���κ�ĳ���
        sin_1=Z_1/L_1;
        cos_1=X_1/L_1;
        rotation_angle=atan2(sin_1,cos_1); % ��Ԫ����ת��
        LI=L_1-Length_1(i); % ��������
        rotation_angle_1=Element_Disp(i,3)-rotation_angle; % ��Ԫi�ڵ�ĵ���ת��
        rotation_angle_2=Element_Disp(i,6)-rotation_angle; % ��Ԫj�ڵ�ĵ���ת��
        Element_Disp(i,:)=[0,0,rotation_angle_1,LI,0,rotation_angle_2]; % ����λ�ƾ���ȥ������λ�ƣ�ʹ������е��Ա���
        KE=beam2d_stiffness_local(E,A,I,Length_1(i),Element_Load_Geo(i,:)); % Length_1��Ӧ�ű���ǰ
        % ���ֲ���������������ò����Ľڵ���
        Element_Load(i,:)=KE*Element_Disp(i,:)';
        Element_Load_Geo(i,:)=Element_Load_Geo(i,:)+Element_Load(i,:);
        % C2״̬
        Coordinate_2(i,:)=Update_Node(N2,:)-Update_Node(N1,:); % ��ԪC1״̬���κ�����
        Length_2(i)=norm(Coordinate_2(i,:));  % ��Ԫ����ǰ�ĳ���
        cs_2=Coordinate_2(i,:)/Length_2(i); % ��Ԫת�ǵ�cos��sinֵ
        T_2=[cs_2(1,1),cs_2(1,2),0,0,0,0;...
            -cs_2(1,2),cs_2(1,1),0,0,0,0;...
            0,0,1,0,0,0;...
            0,0,0,cs_2(1,1),cs_2(1,2),0;...
            0,0,0,-cs_2(1,2),cs_2(1,1),0;...
            0,0,0,0,0,1];
        Element_Load_2(:)=T_2'*Element_Load_Geo(i,:)'; % �õ���Ԫ����������
        Ex_load_global(3*N1-2:3*N1,1)=Element_Load_2(1:3); % i�ڵ����
        Ex_load_global(3*N2-2:3*N2,1)=Element_Load_2(4:6); % j�ڵ����
        Internal_force=Internal_force+Ex_load_global; % �õ�C1״̬�¶�Ӧ��P
    end
    Val=Internal_force-Ex_Load;      % ��������ɺ�����������ؼ�Ĳв�
    Correct_dis=zeros(Num_node,3);   % �����µĽڵ�λ������������ÿ�θ���
    Correct_dis_1=zeros(Num_node,3); % �����µĽڵ�λ�����������ڵ���λ������
    Num_node_dis=size(dis,1);        % δ��Լ���Ľڵ�λ����Ŀ
    k=0;
    while norm(Val)>1e-7 && k<200 || Update_Node(Num_node,1)==0 % ���������������е�����ֹ����
    	    k=k+1;
    	    K=zeros(Num_node*3,Num_node*3);
    	    for i=1:Num_element
    	        N1=Element(i,1);   % ��ȡ��Ԫi�ڵ���
              N2=Element(i,2);   % ��ȡ��Ԫj�ڵ���
              N_E=Element(i,3);  % ��ȡ��Ԫ���
              E=Property(N_E,1); % ��ȡ��Ԫ���ϵ���ģ������λΪN/m^2
              A=Property(N_E,2); % ��ȡ��Ԫ�����������λΪm^2
              I=Property(N_E,3); % ��ȡ��Ԫ���Ծأ���λΪm^4
              Coordinate(i,:)=Update_Node(N2,:)-Update_Node(N1,:); % (a1��)�������ɵĵ�Ԫ��������
              Length1(i)=norm(Coordinate(i,:)); % �������ɵĵ�Ԫ���ȱ仯
              cs=Coordinate(i,:)/Length1(i); % ��Ԫת�ǵ�cos��sinֵ
              KE=beam2d_stiffness(E,A,I,Length1(i),cs,Element_Load_Geo(i,:));
              K=K+Assemble(KE,Element(i,:),Num_node); % �õ��ܸ�
          end
          % �ں����������в���Ϊ����������½���λ�����
          Equation_save=[K,Val];  % �洢��������ϵ��
          Disp_correct=ones(Num_node*3,1);  % ����λ�ƾ������õ�λ������
          Disp_correct(correct,:)=0;        % ����Լ����λ�Ƶ�ֵΪ0�������Ϊ1
          % ����߽�������������
          Equation_save(correct,:)=[];      % ���������ܵ�Լ����λ�ƶ�Ӧ���л�ȥ
          Equation_save(:,correct)=[];      % ���������ܵ�Լ����λ�ƶ�Ӧ���л�ȥ
          Num_column=size(Equation_save,2); % ��ȡ����ϵ���е�����
          % ��ⷽ���飬�õ��ڵ�λ��
          dis_2=-(Equation_save(:,1:Num_column-1))\Equation_save(:,Num_column);
          aa=1;
          for i=1:Num_node*3 % ��λ����ȡ��д��λ���������У��γ��������λ��������
              if Disp_correct(i,1)==1
        	       Disp_correct(i,1)=dis_2(aa,1);
        	       aa=aa+1;
              end	 
          end
          for i=1:Num_node
              Correct_dis(i,:)=Disp_correct(3*i-2:3*i,1); % ���������λ��д��λ�ƾ���
          end
          Correct_dis_1=Correct_dis_1+Correct_dis; % ����λ�ƾ��󣬵õ���λ������
          Internal_force=zeros(Num_node*3,1); % ���ɿյĽڵ������洢����
          Update_Node_1=Update_Node; % C1״̬���ò�����ǰ��λ��
          Update_Node=Node+Displacement(:,1:2)+Correct_dis_1(:,1:2); % C2״̬���ò��������λ��
          % if abs(Update_Node((Num_node-1)/2+1,2))<=1e-7 && abs(Update_Node((Num_node-1)/2+1,1))<=1e-7
          %    break
          % end
          for i=1:Num_element
              Ex_load_global=zeros(Num_node*3,1); % ����ȫ�ֽڵ��������
              for j=1:2
                  % ������λ�ƾ����еõ���ǰ��Ԫ�Ľڵ�λ������
                  Ele_dis(j*3-2:j*3)=Correct_dis(Element(i,j),:);
              end
              N1=Element(i,1);   % ��ȡ��Ԫi�ڵ���
              N2=Element(i,2);   % ��ȡ��Ԫj�ڵ���
              % C1״̬
              Coordinate_1(i,:)=Update_Node_1(N2,:)-Update_Node_1(N1,:); % ��Ԫ����ǰ����
              Length_1(i)=norm(Coordinate_1(i,:));  % ��Ԫ����ǰ�ĳ���
              cs_1=Coordinate_1(i,:)/Length_1(i); % ��Ԫת�ǵ�cos��sinֵ
              T_1=[cs_1(1,1),cs_1(1,2),0,0,0,0;...
                  -cs_1(1,2),cs_1(1,1),0,0,0,0;...
                  0,0,1,0,0,0;...
                  0,0,0,cs_1(1,1),cs_1(1,2),0;...
                  0,0,0,-cs_1(1,2),cs_1(1,1),0;...
                  0,0,0,0,0,1];
              Element_Disp(i,:)=T_1*Ele_dis(:);
              X_1=Length_1(i)+Element_Disp(i,4)-Element_Disp(i,1);
              Z_1=Element_Disp(i,5)-Element_Disp(i,2);
              L_1=(X_1^2+Z_1^2)^0.5; % ��Ԫ���κ�ĳ���
              sin_1=Z_1/L_1;
              cos_1=X_1/L_1;
              rotation_angle=atan2(sin_1,cos_1); % ��Ԫ����ת��
              LI=L_1-Length_1(i); % ��������
              rotation_angle_1=Element_Disp(i,3)-rotation_angle; % ��Ԫi�ڵ�ĵ���ת��
              rotation_angle_2=Element_Disp(i,6)-rotation_angle; % ��Ԫj�ڵ�ĵ���ת��
              Element_Disp(i,:)=[0,0,rotation_angle_1,LI,0,rotation_angle_2]; % ����λ�ƾ���ȥ������λ�ƣ�ʹ������е��Ա���
              KE=beam2d_stiffness_local(E,A,I,Length_1(i),Element_Load_Geo(i,:)); % Length_1��Ӧ�ű���ǰ
              % ���ֲ���������������ò����Ľڵ���
              Element_Load(i,:)=KE*Element_Disp(i,:)';
              Element_Load_Geo(i,:)=Element_Load_Geo(i,:)+Element_Load(i,:);
              % C2״̬
              Coordinate_2(i,:)=Update_Node(N2,:)-Update_Node(N1,:); % ��ԪC1״̬���κ�����
              Length_2(i)=norm(Coordinate_2(i,:));  % ��Ԫ����ǰ�ĳ���
              cs_2=Coordinate_2(i,:)/Length_2(i); % ��Ԫת�ǵ�cos��sinֵ
              T_2=[cs_2(1,1),cs_2(1,2),0,0,0,0;...
                  -cs_2(1,2),cs_2(1,1),0,0,0,0;...
                  0,0,1,0,0,0;...
                  0,0,0,cs_2(1,1),cs_2(1,2),0;...
                  0,0,0,-cs_2(1,2),cs_2(1,1),0;...
                  0,0,0,0,0,1];
              Element_Load_2(:)=T_2'*Element_Load_Geo(i,:)'; % �õ���Ԫ����������
              Ex_load_global(3*N1-2:3*N1,1)=Element_Load_2(1:3); % i�ڵ����
              Ex_load_global(3*N2-2:3*N2,1)=Element_Load_2(4:6); % j�ڵ����
              Internal_force=Internal_force+Ex_load_global; % �õ�C1״̬�¶�Ӧ��P
          end
          Val=Internal_force-Ex_Load;      % ��������ɺ�����������ؼ�Ĳв�
          Val(correct,:)=0;
    end
    Displacement=Displacement+Correct_dis_1;
    x(n+1,1)=Displacement(Num_node,2);
    y(n+1,1)=Ex_Load(Num_node*3,1);
end

figure
plot(x,y,'r','LineWidth',2)
title('���������ص����λ������')
xlabel('λ�ƣ�m��')
ylabel('���أ�N.m��')
grid on

figure
plot(Update_Node(1:Num_node,1),Update_Node(1:Num_node,2))