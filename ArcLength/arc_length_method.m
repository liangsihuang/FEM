% 弧长法(带有超平面更新的超平面约束) 主程序 
Node=xlsread('Node'); % 输入节点坐标信息
Element=xlsread('Element'); % 输入单元编号及单元两侧节点编号
Boundary=xlsread('Boundary'); % 输入边界条件
Section=xlsread('Section'); % 输入截面参数，第一列为弹性模量，第二列为截面面积，第三列为截面惯性矩
Force=xlsread('Force'); % 输入外荷载，形成荷载信息矩阵
Huchang=0.05; % 初始弧长
N_Node=size(Node,1); % 节点数量
N_Element=size(Element,1); % 单元数量
N_Force=size(Force,1);% 外部作用荷载数量
N_Boundary=size(Boundary,1); % 约束数量
Update_Node=Node; % 更新的节点坐标,x坐标与y坐标
All_Disp=zeros(N_Node,3); % 生成位移存储矩阵，初始位移均为0
All_F=zeros(N_Node*3,1);  % 总荷载向量
zlamda=0; % 总荷载系数
Element_F=zeros(N_Element,6); % 单元节点外荷载
Force_Transfer=zeros(N_Node*3,1); % 外荷载矩阵
Internal_F=zeros(N_Node*3,1); % 空的节点内力存储矩阵
C=zeros(N_Element,2); % 单元节点坐标增量存储矩阵
L=zeros(N_Element,1); % 单元长度存储矩阵
Y=zeros(1,1); % 绘图变量
X=zeros(1,1); % 绘图变量
for i=1:N_Force
    Force_Transfer(Force(i,1)*3-2:Force(i,1)*3,1)=Force(i,2:4)'; % 将荷载信息矩阵转换为单列的外荷载矩阵
end
z_lamda=0;
n=0; % 增量步步数记录
while All_Disp(6,2)>-2 % 弧长法求解终止条件：拱跨中位移大于限值
    n=n+1;
    K_Global=zeros(N_Node*3,N_Node*3); % 形成空的总刚矩阵
    for i=1:N_Element    % 生成总刚
        N1=Element(i,1); % 读取单元i节点编号
        N2=Element(i,2); % 读取单元j节点编号
        C(i,:)=Update_Node(N2,:)-Update_Node(N1,:); % 因变形造成的单元坐标增量
        L(i)=norm(C(i,:)); % 因变形造成的单元长度变化
        cs=C(i,:)/L(i);    % 单元转角的cos及sin值
        lamda=lamdabianhuan(cs); % 整体坐标与局部坐标转换矩阵
        K_Element_jb=beam2d_stiffness(i,L(i),Element_F(i,:),Element,Section); % 得到单元在局部坐标下的刚度矩阵
        K_Element_zt=lamda'*K_Element_jb*lamda; % 得到单元在整体坐标下的刚度矩阵
        K_Global=K_Global+Assemble(K_Element_zt,Element(i,:),N_Node); % 组装单刚得到总刚
    end
    [XZK,XZP]=bianjiexiuzheng(Boundary,Force_Transfer,K_Global,N_Node); % 进行边界修正
    dis1=XZK\XZP; % 求解位移
    lamda=Huchang/(dis1'*dis1+1)^0.5; % 得到本增量步第一次迭代后的荷载系数
    Disp_Transfer=lamda*dis1; % 本增量步第一次迭代后的位移 ?????
    r_m=[Disp_Transfer',lamda]'; % 弧长记录
    % 荷载系数符号判定，锐角为正，钝角为负
    if n>1
        if r_m'*R_m<=0
            lamda=-lamda;
            Disp_Transfer=-Disp_Transfer;
            r_m=-r_m;
        end
    end
   R_m=zeros(3*N_Node+1,1); % 超平面弧长储存
   R_m=R_m+r_m; % 超平面弧长的更新
   zlamda=zlamda+lamda; % 总荷载系数更新
   All_F=zlamda*Force_Transfer; % 总荷载更新
   for i=1:N_Node
       All_Disp(i,:)=All_Disp(i,:)+Disp_Transfer(3*i-2:3*i)'; % 总位移更新
   end
   for i=1:N_Element;
       Disp_Element_zt=zeros(6,1); % 整体坐标下单元位移存储矩阵
       N1=Element(i,1); % 读取单元i节点编号
       N2=Element(i,2); % 读取单元j节点编号
       for j=1:3 % 读取单元i节点位移
           Disp_Element_zt(j,1)=Disp_Transfer(3*N1+j-3);
       end
       for j=4:6 % 读取单元i节点位移
           Disp_Element_zt(j,1)=Disp_Transfer(3*N2+j-6);
       end
       C(i,:)=Update_Node(N2,:)-Update_Node(N1,:); % 因变形造成的单元坐标增量
       L(i)=norm(C(i,:)); % 因变形造成的单元长度变化
       cs=C(i,:)/L(i); % 单元转角的cos及sin值
       lamda=lamdabianhuan(cs); % 整体坐标与局部坐标转换矩阵
       Disp_Element_jb=lamda*Disp_Element_zt; % 得到单元在局部坐标下的变形
       Disp_Element_jb_tanxing=tanxingbianxing(Disp_Element_jb,L(i)); % 得到单元在局部坐标下的弹性变形
       K_Element_jb=beam2d_stiffness(i,L(i),Element_F(i,:),Element,Section); % 得到单元在局部坐标下的刚度矩阵
       Element_F_jb=K_Element_jb*Disp_Element_jb_tanxing; % 得到单元在局部坐标下的内力
       Element_F(i,:)=Element_F(i,:)+Element_F_jb(:,1)';  % 更新单元内力
   end
   Update_Node=Node+All_Disp(:,1:2); % 根据位移信息更新的节点坐标,x坐标与y坐标
   Internal_F=zeros(N_Node*3,1); % 空的节点内力存储矩阵
   for i=1:N_Element
       F_Global=zeros(N_Node*3,1); % 生成整体坐标下节点荷载向量
       N1=Element(i,1); % 读取单元i节点编号
       N2=Element(i,2); % 读取单元j节点编号
       C(i,:)=Update_Node(N2,:)-Update_Node(N1,:); % 因变形造成的单元坐标增量
       L(i)=norm(C(i,:)); % 因变形造成的单元长度变化
       cs=C(i,:)/L(i); % 单元转角的cos及sin值
       lamda=lamdabianhuan(cs); % 整体坐标与局部坐标转换矩阵
       Element_All_F_jb=zeros(6,1); % 生成局部节点荷载向量存储矩阵
       Element_All_F_jb(:,1)=Element_F(i,:); % 生成局部节点荷载向量
       Element_All_F_zt=lamda'*Element_All_F_jb; % 得到整体坐标系下节点荷载
       F_Global(3*N1-2:3*N1,1)=Element_All_F_zt(1:3,1); % 得到整体坐标下节点荷载
       F_Global(3*N2-2:3*N2,1)=Element_All_F_zt(4:6,1); % 得到整体坐标下节点荷载
       Internal_F=Internal_F+F_Global; % 节点内力更新
   end
   k=1; % 增量步内迭代次数标记
   eval(['Val',num2str(n),'_',num2str(k),'=All_F-Internal_F;']) % 得到第n个增量步第k次迭代的荷载残差
   for i=1:N_Node % 对荷载残差进行边界修正
       for j=2:4
           if Boundary(i,j)==0
              eval(['Val',num2str(n),'_',num2str(k),'(3*i+j-4)=0;'])
           end
       end
   end
   while eval(['norm(Val',num2str(n),'_',num2str(k),')>1e-5 && k<20']) % 第n个增量步迭代终止条件
        k=k+1;
        K_Global=zeros(N_Node*3,N_Node*3); % 形成空的总刚矩阵
        for i=1:N_Element    % 生成总刚
            N1=Element(i,1); % 读取单元i节点编号
            N2=Element(i,2); % 读取单元j节点编号
            C(i,:)=Update_Node(N2,:)-Update_Node(N1,:); % 因变形造成的单元坐标增量
            L(i)=norm(C(i,:)); % 因变形造成的单元长度变化
            cs=C(i,:)/L(i);    % 单元转角的cos及sin值
            lamda=lamdabianhuan(cs); % 整体坐标与局部坐标转换矩阵
            K_Element_jb=beam2d_stiffness(i,L(i),Element_F(i,:),Element,Section); % 得到单元在局部坐标下的刚度矩阵
            K_Element_zt=lamda'*K_Element_jb*lamda; % 得到单元在整体坐标下的刚度矩阵
            K_Global=K_Global+Assemble(K_Element_zt,Element(i,:),N_Node); % 组装单刚得到总刚
        end
        eval(['[XZK,XZP]=bianjiexiuzheng(Boundary,Val',num2str(n),'_',num2str(k-1),',K_Global,N_Node);']) % 在荷载残差条件下的边界修正
        dis2=XZK\XZP; % 求解荷载残差导致的位移
        [XZK,XZP]=bianjiexiuzheng(Boundary,Force_Transfer,K_Global,N_Node); % 进行边界修正
        dis3=XZK\XZP; % 求解外荷载导致的位移
        lamda=-dis1'*dis2/(dis1'*dis3+1); % 求解荷载系数
        Disp_Transfer=lamda*dis3+dis2; % 迭代后的位移
        r_k=[Disp_Transfer',lamda]'; % 弧长记录
        R_m=R_m+r_k; % 弧长更新
        zlamda=zlamda+lamda; % 总荷载系数更新
        All_F=zlamda*Force_Transfer; % 总荷载更新
        for i=1:N_Node
            All_Disp(i,:)=All_Disp(i,:)+Disp_Transfer(3*i-2:3*i)'; % 总位移更新
        end
        for i=1:N_Element;
            Disp_Element_zt=zeros(6,1); % 整体坐标下单元位移存储矩阵
            N1=Element(i,1); % 读取单元i节点编号
            N2=Element(i,2); % 读取单元j节点编号
            for j=1:3 % 读取单元i节点位移
                Disp_Element_zt(j,1)=Disp_Transfer(3*N1+j-3);
            end
            for j=4:6 % 读取单元i节点位移
                Disp_Element_zt(j,1)=Disp_Transfer(3*N2+j-6);
            end
            C(i,:)=Update_Node(N2,:)-Update_Node(N1,:); % 因变形造成的单元坐标增量
            L(i)=norm(C(i,:)); % 因变形造成的单元长度变化
            cs=C(i,:)/L(i); % 单元转角的cos及sin值
            lamda=lamdabianhuan(cs); % 整体坐标与局部坐标转换矩阵
            Disp_Element_jb=lamda*Disp_Element_zt; % 得到单元在局部坐标下的变形
            Disp_Element_jb_tanxing=tanxingbianxing(Disp_Element_jb,L(i)); % 得到单元在局部坐标下的弹性变形
            K_Element_jb=beam2d_stiffness(i,L(i),Element_F(i,:),Element,Section); % 得到单元在局部坐标下的刚度矩阵
            Element_F_jb=K_Element_jb*Disp_Element_jb_tanxing; % 得到单元在局部坐标下的内力
            Element_F(i,:)=Element_F(i,:)+Element_F_jb(:,1)';  % 更新单元内力
        end
        Update_Node=Node+All_Disp(:,1:2); % 根据位移信息更新的节点坐标,x坐标与y坐标
        Internal_F=zeros(N_Node*3,1); % 空的节点内力存储矩阵
        for i=1:N_Element
            F_Global=zeros(N_Node*3,1); % 生成整体坐标下节点荷载向量
            N1=Element(i,1); % 读取单元i节点编号
            N2=Element(i,2); % 读取单元j节点编号
            C(i,:)=Update_Node(N2,:)-Update_Node(N1,:); % 因变形造成的单元坐标增量
            L(i)=norm(C(i,:)); % 因变形造成的单元长度变化
            cs=C(i,:)/L(i); % 单元转角的cos及sin值
            lamda=lamdabianhuan(cs); % 整体坐标与局部坐标转换矩阵
            Element_All_F_jb=zeros(6,1); % 生成局部节点荷载向量存储矩阵
            Element_All_F_jb(:,1)=Element_F(i,:); % 生成局部节点荷载向量
            Element_All_F_zt=lamda'*Element_All_F_jb; % 得到整体坐标系下节点荷载
            F_Global(3*N1-2:3*N1,1)=Element_All_F_zt(1:3,1); % 得到整体坐标下节点荷载
            F_Global(3*N2-2:3*N2,1)=Element_All_F_zt(4:6,1); % 得到整体坐标下节点荷载
            Internal_F=Internal_F+F_Global; % 节点内力更新
        end
        eval(['Val',num2str(n),'_',num2str(k),'=All_F-Internal_F;']) % 得到第n个增量步第k次迭代的荷载残差
        for i=1:N_Node % 对荷载残差进行边界修正
            for j=2:4
                if Boundary(i,j)==0
                   eval(['Val',num2str(n),'_',num2str(k),'(3*i+j-4)=0;'])
                end
            end
        end
   end
   Y=[Y,0];
   X=[X,0];
   Y(1,n+1)=-All_F(17,1);   % 记录跨中荷载
   X(1,n+1)=-All_Disp(6,2); % 记录跨中竖向位移
end
% 绘图
plot(X,Y);
grid on;