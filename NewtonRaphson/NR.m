% Newton_Raphson法 主程序 
clear;
clc;
% 参数输入
Node=xlsread('Data.xlsx','Node'); % 输入节点坐标信息
Element=xlsread('Data.xlsx','Element'); % 输入单元编号及单元两侧节点编号
Boundary=xlsread('Data.xlsx','Boundary'); % 输入边界条件（本题中仅1节点固端约束）
Property=xlsread('Data.xlsx','Property'); % 输入参数，第一列为弹性模量，第二列为截面面积，第三列为截面惯性矩
External_Load=xlsread('Data.xlsx','External_Load'); % 输入外荷载，形成荷载信息矩阵（本题仅有作用于悬臂梁端部集中弯矩）

Num_node=size(Node,1); % 节点数量
Num_element=size(Element,1); % 单元数量
Num_el=size(External_Load,1); % 外部作用荷载数量
Num_bn=size(Boundary,1); % 有约束节点数量
Displacement=zeros(Num_node,3); % 生成位移存储矩阵，初始位移均为0
Internal_force=zeros(Num_node*3,1); % 生成空的节点内力存储矩阵
Dis=zeros(Num_node,3); % 生成空的节点位移转角迭代存储矩阵
Ex_Load=zeros(Num_node*3,1); % 生成空的初始节点外荷载存储矩阵
Total_load=zeros(Num_node*3,1); % 总荷载向量
Total_step=1000; % 迭代总次数

for i=1:Num_el
    Total_load(External_Load(i,1)*3-2:External_Load(i,1)*3,1)=External_Load(i,2:4)'; % 将荷载信息矩阵转换为单列的外荷载矩阵
end

% 边界条件修正信息，求出约束节点的标号，便于刚度、荷载矩阵归0
correct=[]; % 生成空矩阵
for i=1:Num_bn
    if Boundary(i,2)==0
       m=Boundary(i,1)*3-2;    % 受约束节点x向位移为0时所对应的行列编号为 3*节点编号-2
       correct=[correct,[m]];  % 记录行列编号
    end
    if Boundary(i,3)==0
       m=Boundary(i,1)*3-1;    % 受约束节点y向位移为0时所对应的行列编号为 3*节点编号-1
       correct=[correct,[m]];  % 记录行列编号 
    end
    if Boundary(i,4)==0
       m=Boundary(i,1)*3;      % 受约束节点转角为0时所对应的行列编号为 3*节点编号
       correct=[correct,[m]];  % 记录行列编号
    end
end

Update_Node=Node+Dis(:,1:2); % 根据位移信息更新的节点坐标,x坐标与y坐标
Coordinate=zeros(Num_node,2); % 单元节点坐标增量存储矩阵
Length1=zeros(Num_node,1); % 单元长度储存
Element_Load=zeros(Num_element,6); % 单元节点外荷载
Element_Disp=zeros(Num_element,6); % 单元节点位移
Element_Load_Geo=zeros(Num_element,6); % 单元节点外荷载中与几何矩阵有关的外荷载
Element_Load_2=zeros(1,6); % 单个单元节点力临时存储矩阵
Ele_dis=zeros(1,6); % 单个单元节点位移临时存储矩阵
x(1,1)=0; % 绘图变量
y(1,1)=0; % 绘图变量

for n=1:Total_step
    n=n+1;
    K=zeros(Num_node*3,Num_node*3); % 形成空的总刚矩阵
    for i=1:Num_element    % 生成总刚 
        N1=Element(i,1);   % 读取单元i节点编号
        N2=Element(i,2);   % 读取单元j节点编号
        N_E=Element(i,3);  % 读取单元编号
        E=Property(N_E,1); % 读取单元材料弹性模量，单位为N/m^2
        A=Property(N_E,2); % 读取单元截面面积，单位为m^2
        I=Property(N_E,3); % 读取单元惯性矩，单位为m^4
        Coordinate(i,:)=Update_Node(N2,:)-Update_Node(N1,:); % 因变形造成的单元坐标增量
        Length1(i)=norm(Coordinate(i,:)); % 因变形造成的单元长度变化
        cs=Coordinate(i,:)/Length1(i); % 单元转角的cos及sin值
        KE=beam2d_stiffness(E,A,I,Length1(i),cs,Element_Load_Geo(i,:));
        K=K+Assemble(KE,Element(i,:),Num_node); % 得到总刚
    end
    Load_step=Total_load/Total_step;  % 计算荷载初始增量步大小
    Equation_save=[K,Load_step];      % 存储增量方程系数
    Disp_correct=ones(Num_node*3,1);  % 修正位移矩阵所用的位移向量
    Disp_correct(correct,:)=0;        % 令受约束的位移的值为0，其余的为1
    % 引入边界条件修正方程
    Equation_save(correct,:)=[];      % 将方程中受到约束的位移对应的行划去
    Equation_save(:,correct)=[];      % 将方程中受到约束的位移对应的列划去
    Num_column=size(Equation_save,2); % 读取方程系数中的列数
    % 求解方程组，得到节点位移
    dis=Equation_save(:,1:Num_column-1)\Equation_save(:,Num_column);
    aa=1;
    for i=1:Num_node*3 % 将位移提取并写入位移列向量中，形成修正后的位移列向量
        if Disp_correct(i,1)==1
        	 Disp_correct(i,1)=dis(aa,1);
        	 aa=aa+1;
        end	 
    end
    for i=1:Num_node
        Dis(i,:)=Disp_correct(3*i-2:3*i,1); % 将修正后的位移写入位移矩阵
    end
    Displacement=Dis+Displacement; % 更新位移矩阵，得到总位移增量
    Ex_Load=Ex_Load+Load_step; % 更新外荷载向量
    Internal_force=zeros(Num_node*3,1); % 生成空的节点内力存储矩阵
    Update_Node_1=Update_Node; % C1状态，该步迭代前的位置
    Update_Node=Node+Displacement(:,1:2); % C2状态，该步迭代后的位置
    for i=1:Num_element
        Ex_load_global=zeros(Num_node*3,1); % 生成全局节点荷载向量
        for j=1:2
            % 从整体位移矩阵中得到当前单元的节点位移向量
            Ele_dis(j*3-2:j*3)=Dis(Element(i,j),:);
        end
        N1=Element(i,1);   % 读取单元i节点编号
        N2=Element(i,2);   % 读取单元j节点编号
        % C1状态
        Coordinate_1(i,:)=Update_Node_1(N2,:)-Update_Node_1(N1,:); % 单元变形前坐标
        Length_1(i)=norm(Coordinate_1(i,:));  % 单元变形前的长度
        cs_1=Coordinate_1(i,:)/Length_1(i); % 单元转角的cos及sin值
        T_1=[cs_1(1,1),cs_1(1,2),0,0,0,0;...
            -cs_1(1,2),cs_1(1,1),0,0,0,0;...
            0,0,1,0,0,0;...
            0,0,0,cs_1(1,1),cs_1(1,2),0;...
            0,0,0,-cs_1(1,2),cs_1(1,1),0;...
            0,0,0,0,0,1];
        Element_Disp(i,:)=T_1*Ele_dis(:);
        X_1=Length_1(i)+Element_Disp(i,4)-Element_Disp(i,1);
        Z_1=Element_Disp(i,5)-Element_Disp(i,2);
        L_1=(X_1^2+Z_1^2)^0.5; % 单元变形后的长度
        sin_1=Z_1/L_1;
        cos_1=X_1/L_1;
        rotation_angle=atan2(sin_1,cos_1); % 单元刚体转角
        LI=L_1-Length_1(i); % 长度增量
        rotation_angle_1=Element_Disp(i,3)-rotation_angle; % 单元i节点的弹性转角
        rotation_angle_2=Element_Disp(i,6)-rotation_angle; % 单元j节点的弹性转角
        Element_Disp(i,:)=[0,0,rotation_angle_1,LI,0,rotation_angle_2]; % 更新位移矩阵，去除刚体位移，使其仅含有弹性变形
        KE=beam2d_stiffness_local(E,A,I,Length_1(i),Element_Load_Geo(i,:)); % Length_1对应着变形前
        % 求解局部坐标下因荷载作用产生的节点力
        Element_Load(i,:)=KE*Element_Disp(i,:)';
        Element_Load_Geo(i,:)=Element_Load_Geo(i,:)+Element_Load(i,:);
        % C2状态
        Coordinate_2(i,:)=Update_Node(N2,:)-Update_Node(N1,:); % 单元C1状态变形后坐标
        Length_2(i)=norm(Coordinate_2(i,:));  % 单元变形前的长度
        cs_2=Coordinate_2(i,:)/Length_2(i); % 单元转角的cos及sin值
        T_2=[cs_2(1,1),cs_2(1,2),0,0,0,0;...
            -cs_2(1,2),cs_2(1,1),0,0,0,0;...
            0,0,1,0,0,0;...
            0,0,0,cs_2(1,1),cs_2(1,2),0;...
            0,0,0,-cs_2(1,2),cs_2(1,1),0;...
            0,0,0,0,0,1];
        Element_Load_2(:)=T_2'*Element_Load_Geo(i,:)'; % 得到单元荷载行向量
        Ex_load_global(3*N1-2:3*N1,1)=Element_Load_2(1:3); % i节点荷载
        Ex_load_global(3*N2-2:3*N2,1)=Element_Load_2(4:6); % j节点荷载
        Internal_force=Internal_force+Ex_load_global; % 得到C1状态下对应的P
    end
    Val=Internal_force-Ex_Load;      % 求解迭代完成后内力与外荷载间的残差
    Correct_dis=zeros(Num_node,3);   % 构造新的节点位移向量，用于每次更新
    Correct_dis_1=zeros(Num_node,3); % 构造新的节点位移向量，用于叠加位移增量
    Num_node_dis=size(dis,1);        % 未受约束的节点位移数目
    k=0;
    while norm(Val)>1e-7 && k<200 || Update_Node(Num_node,1)==0 % 单个荷载增量步中迭代终止条件
    	    k=k+1;
    	    K=zeros(Num_node*3,Num_node*3);
    	    for i=1:Num_element
    	        N1=Element(i,1);   % 读取单元i节点编号
              N2=Element(i,2);   % 读取单元j节点编号
              N_E=Element(i,3);  % 读取单元编号
              E=Property(N_E,1); % 读取单元材料弹性模量，单位为N/m^2
              A=Property(N_E,2); % 读取单元截面面积，单位为m^2
              I=Property(N_E,3); % 读取单元惯性矩，单位为m^4
              Coordinate(i,:)=Update_Node(N2,:)-Update_Node(N1,:); % (a1下)因变形造成的单元坐标增量
              Length1(i)=norm(Coordinate(i,:)); % 因变形造成的单元长度变化
              cs=Coordinate(i,:)/Length1(i); % 单元转角的cos及sin值
              KE=beam2d_stiffness(E,A,I,Length1(i),cs,Element_Load_Geo(i,:));
              K=K+Assemble(KE,Element(i,:),Num_node); % 得到总刚
          end
          % 在荷载与内力残差作为外荷载条件下进行位移求解
          Equation_save=[K,Val];  % 存储增量方程系数
          Disp_correct=ones(Num_node*3,1);  % 修正位移矩阵所用的位移向量
          Disp_correct(correct,:)=0;        % 令受约束的位移的值为0，其余的为1
          % 引入边界条件修正方程
          Equation_save(correct,:)=[];      % 将方程中受到约束的位移对应的行划去
          Equation_save(:,correct)=[];      % 将方程中受到约束的位移对应的列划去
          Num_column=size(Equation_save,2); % 读取方程系数中的列数
          % 求解方程组，得到节点位移
          dis_2=-(Equation_save(:,1:Num_column-1))\Equation_save(:,Num_column);
          aa=1;
          for i=1:Num_node*3 % 将位移提取并写入位移列向量中，形成修正后的位移列向量
              if Disp_correct(i,1)==1
        	       Disp_correct(i,1)=dis_2(aa,1);
        	       aa=aa+1;
              end	 
          end
          for i=1:Num_node
              Correct_dis(i,:)=Disp_correct(3*i-2:3*i,1); % 将修正后的位移写入位移矩阵
          end
          Correct_dis_1=Correct_dis_1+Correct_dis; % 更新位移矩阵，得到总位移增量
          Internal_force=zeros(Num_node*3,1); % 生成空的节点内力存储矩阵
          Update_Node_1=Update_Node; % C1状态，该步迭代前的位置
          Update_Node=Node+Displacement(:,1:2)+Correct_dis_1(:,1:2); % C2状态，该步迭代后的位置
          % if abs(Update_Node((Num_node-1)/2+1,2))<=1e-7 && abs(Update_Node((Num_node-1)/2+1,1))<=1e-7
          %    break
          % end
          for i=1:Num_element
              Ex_load_global=zeros(Num_node*3,1); % 生成全局节点荷载向量
              for j=1:2
                  % 从整体位移矩阵中得到当前单元的节点位移向量
                  Ele_dis(j*3-2:j*3)=Correct_dis(Element(i,j),:);
              end
              N1=Element(i,1);   % 读取单元i节点编号
              N2=Element(i,2);   % 读取单元j节点编号
              % C1状态
              Coordinate_1(i,:)=Update_Node_1(N2,:)-Update_Node_1(N1,:); % 单元变形前坐标
              Length_1(i)=norm(Coordinate_1(i,:));  % 单元变形前的长度
              cs_1=Coordinate_1(i,:)/Length_1(i); % 单元转角的cos及sin值
              T_1=[cs_1(1,1),cs_1(1,2),0,0,0,0;...
                  -cs_1(1,2),cs_1(1,1),0,0,0,0;...
                  0,0,1,0,0,0;...
                  0,0,0,cs_1(1,1),cs_1(1,2),0;...
                  0,0,0,-cs_1(1,2),cs_1(1,1),0;...
                  0,0,0,0,0,1];
              Element_Disp(i,:)=T_1*Ele_dis(:);
              X_1=Length_1(i)+Element_Disp(i,4)-Element_Disp(i,1);
              Z_1=Element_Disp(i,5)-Element_Disp(i,2);
              L_1=(X_1^2+Z_1^2)^0.5; % 单元变形后的长度
              sin_1=Z_1/L_1;
              cos_1=X_1/L_1;
              rotation_angle=atan2(sin_1,cos_1); % 单元刚体转角
              LI=L_1-Length_1(i); % 长度增量
              rotation_angle_1=Element_Disp(i,3)-rotation_angle; % 单元i节点的弹性转角
              rotation_angle_2=Element_Disp(i,6)-rotation_angle; % 单元j节点的弹性转角
              Element_Disp(i,:)=[0,0,rotation_angle_1,LI,0,rotation_angle_2]; % 更新位移矩阵，去除刚体位移，使其仅含有弹性变形
              KE=beam2d_stiffness_local(E,A,I,Length_1(i),Element_Load_Geo(i,:)); % Length_1对应着变形前
              % 求解局部坐标下因荷载作用产生的节点力
              Element_Load(i,:)=KE*Element_Disp(i,:)';
              Element_Load_Geo(i,:)=Element_Load_Geo(i,:)+Element_Load(i,:);
              % C2状态
              Coordinate_2(i,:)=Update_Node(N2,:)-Update_Node(N1,:); % 单元C1状态变形后坐标
              Length_2(i)=norm(Coordinate_2(i,:));  % 单元变形前的长度
              cs_2=Coordinate_2(i,:)/Length_2(i); % 单元转角的cos及sin值
              T_2=[cs_2(1,1),cs_2(1,2),0,0,0,0;...
                  -cs_2(1,2),cs_2(1,1),0,0,0,0;...
                  0,0,1,0,0,0;...
                  0,0,0,cs_2(1,1),cs_2(1,2),0;...
                  0,0,0,-cs_2(1,2),cs_2(1,1),0;...
                  0,0,0,0,0,1];
              Element_Load_2(:)=T_2'*Element_Load_Geo(i,:)'; % 得到单元荷载行向量
              Ex_load_global(3*N1-2:3*N1,1)=Element_Load_2(1:3); % i节点荷载
              Ex_load_global(3*N2-2:3*N2,1)=Element_Load_2(4:6); % j节点荷载
              Internal_force=Internal_force+Ex_load_global; % 得到C1状态下对应的P
          end
          Val=Internal_force-Ex_Load;      % 求解迭代完成后内力与外荷载间的残差
          Val(correct,:)=0;
    end
    Displacement=Displacement+Correct_dis_1;
    x(n+1,1)=Displacement(Num_node,2);
    y(n+1,1)=Ex_Load(Num_node*3,1);
end

figure
plot(x,y,'r','LineWidth',2)
title('悬臂梁加载点荷载位移曲线')
xlabel('位移（m）')
ylabel('荷载（N.m）')
grid on

figure
plot(Update_Node(1:Num_node,1),Update_Node(1:Num_node,2))