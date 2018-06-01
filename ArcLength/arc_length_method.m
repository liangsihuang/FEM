% ������(���г�ƽ����µĳ�ƽ��Լ��) ������ 
Node=xlsread('Node'); % ����ڵ�������Ϣ
Element=xlsread('Element'); % ���뵥Ԫ��ż���Ԫ����ڵ���
Boundary=xlsread('Boundary'); % ����߽�����
Section=xlsread('Section'); % ��������������һ��Ϊ����ģ�����ڶ���Ϊ���������������Ϊ������Ծ�
Force=xlsread('Force'); % ��������أ��γɺ�����Ϣ����
Huchang=0.05; % ��ʼ����
N_Node=size(Node,1); % �ڵ�����
N_Element=size(Element,1); % ��Ԫ����
N_Force=size(Force,1);% �ⲿ���ú�������
N_Boundary=size(Boundary,1); % Լ������
Update_Node=Node; % ���µĽڵ�����,x������y����
All_Disp=zeros(N_Node,3); % ����λ�ƴ洢���󣬳�ʼλ�ƾ�Ϊ0
All_F=zeros(N_Node*3,1);  % �ܺ�������
zlamda=0; % �ܺ���ϵ��
Element_F=zeros(N_Element,6); % ��Ԫ�ڵ������
Force_Transfer=zeros(N_Node*3,1); % ����ؾ���
Internal_F=zeros(N_Node*3,1); % �յĽڵ������洢����
C=zeros(N_Element,2); % ��Ԫ�ڵ����������洢����
L=zeros(N_Element,1); % ��Ԫ���ȴ洢����
Y=zeros(1,1); % ��ͼ����
X=zeros(1,1); % ��ͼ����
for i=1:N_Force
    Force_Transfer(Force(i,1)*3-2:Force(i,1)*3,1)=Force(i,2:4)'; % ��������Ϣ����ת��Ϊ���е�����ؾ���
end
z_lamda=0;
n=0; % ������������¼
while All_Disp(6,2)>-2 % �����������ֹ������������λ�ƴ�����ֵ
    n=n+1;
    K_Global=zeros(N_Node*3,N_Node*3); % �γɿյ��ܸվ���
    for i=1:N_Element    % �����ܸ�
        N1=Element(i,1); % ��ȡ��Ԫi�ڵ���
        N2=Element(i,2); % ��ȡ��Ԫj�ڵ���
        C(i,:)=Update_Node(N2,:)-Update_Node(N1,:); % �������ɵĵ�Ԫ��������
        L(i)=norm(C(i,:)); % �������ɵĵ�Ԫ���ȱ仯
        cs=C(i,:)/L(i);    % ��Ԫת�ǵ�cos��sinֵ
        lamda=lamdabianhuan(cs); % ����������ֲ�����ת������
        K_Element_jb=beam2d_stiffness(i,L(i),Element_F(i,:),Element,Section); % �õ���Ԫ�ھֲ������µĸնȾ���
        K_Element_zt=lamda'*K_Element_jb*lamda; % �õ���Ԫ�����������µĸնȾ���
        K_Global=K_Global+Assemble(K_Element_zt,Element(i,:),N_Node); % ��װ���յõ��ܸ�
    end
    [XZK,XZP]=bianjiexiuzheng(Boundary,Force_Transfer,K_Global,N_Node); % ���б߽�����
    dis1=XZK\XZP; % ���λ��
    lamda=Huchang/(dis1'*dis1+1)^0.5; % �õ�����������һ�ε�����ĺ���ϵ��
    Disp_Transfer=lamda*dis1; % ����������һ�ε������λ�� ?????
    r_m=[Disp_Transfer',lamda]'; % ������¼
    % ����ϵ�������ж������Ϊ�����۽�Ϊ��
    if n>1
        if r_m'*R_m<=0
            lamda=-lamda;
            Disp_Transfer=-Disp_Transfer;
            r_m=-r_m;
        end
    end
   R_m=zeros(3*N_Node+1,1); % ��ƽ�满������
   R_m=R_m+r_m; % ��ƽ�满���ĸ���
   zlamda=zlamda+lamda; % �ܺ���ϵ������
   All_F=zlamda*Force_Transfer; % �ܺ��ظ���
   for i=1:N_Node
       All_Disp(i,:)=All_Disp(i,:)+Disp_Transfer(3*i-2:3*i)'; % ��λ�Ƹ���
   end
   for i=1:N_Element;
       Disp_Element_zt=zeros(6,1); % ���������µ�Ԫλ�ƴ洢����
       N1=Element(i,1); % ��ȡ��Ԫi�ڵ���
       N2=Element(i,2); % ��ȡ��Ԫj�ڵ���
       for j=1:3 % ��ȡ��Ԫi�ڵ�λ��
           Disp_Element_zt(j,1)=Disp_Transfer(3*N1+j-3);
       end
       for j=4:6 % ��ȡ��Ԫi�ڵ�λ��
           Disp_Element_zt(j,1)=Disp_Transfer(3*N2+j-6);
       end
       C(i,:)=Update_Node(N2,:)-Update_Node(N1,:); % �������ɵĵ�Ԫ��������
       L(i)=norm(C(i,:)); % �������ɵĵ�Ԫ���ȱ仯
       cs=C(i,:)/L(i); % ��Ԫת�ǵ�cos��sinֵ
       lamda=lamdabianhuan(cs); % ����������ֲ�����ת������
       Disp_Element_jb=lamda*Disp_Element_zt; % �õ���Ԫ�ھֲ������µı���
       Disp_Element_jb_tanxing=tanxingbianxing(Disp_Element_jb,L(i)); % �õ���Ԫ�ھֲ������µĵ��Ա���
       K_Element_jb=beam2d_stiffness(i,L(i),Element_F(i,:),Element,Section); % �õ���Ԫ�ھֲ������µĸնȾ���
       Element_F_jb=K_Element_jb*Disp_Element_jb_tanxing; % �õ���Ԫ�ھֲ������µ�����
       Element_F(i,:)=Element_F(i,:)+Element_F_jb(:,1)';  % ���µ�Ԫ����
   end
   Update_Node=Node+All_Disp(:,1:2); % ����λ����Ϣ���µĽڵ�����,x������y����
   Internal_F=zeros(N_Node*3,1); % �յĽڵ������洢����
   for i=1:N_Element
       F_Global=zeros(N_Node*3,1); % �������������½ڵ��������
       N1=Element(i,1); % ��ȡ��Ԫi�ڵ���
       N2=Element(i,2); % ��ȡ��Ԫj�ڵ���
       C(i,:)=Update_Node(N2,:)-Update_Node(N1,:); % �������ɵĵ�Ԫ��������
       L(i)=norm(C(i,:)); % �������ɵĵ�Ԫ���ȱ仯
       cs=C(i,:)/L(i); % ��Ԫת�ǵ�cos��sinֵ
       lamda=lamdabianhuan(cs); % ����������ֲ�����ת������
       Element_All_F_jb=zeros(6,1); % ���ɾֲ��ڵ���������洢����
       Element_All_F_jb(:,1)=Element_F(i,:); % ���ɾֲ��ڵ��������
       Element_All_F_zt=lamda'*Element_All_F_jb; % �õ���������ϵ�½ڵ����
       F_Global(3*N1-2:3*N1,1)=Element_All_F_zt(1:3,1); % �õ����������½ڵ����
       F_Global(3*N2-2:3*N2,1)=Element_All_F_zt(4:6,1); % �õ����������½ڵ����
       Internal_F=Internal_F+F_Global; % �ڵ���������
   end
   k=1; % �������ڵ����������
   eval(['Val',num2str(n),'_',num2str(k),'=All_F-Internal_F;']) % �õ���n����������k�ε����ĺ��زв�
   for i=1:N_Node % �Ժ��زв���б߽�����
       for j=2:4
           if Boundary(i,j)==0
              eval(['Val',num2str(n),'_',num2str(k),'(3*i+j-4)=0;'])
           end
       end
   end
   while eval(['norm(Val',num2str(n),'_',num2str(k),')>1e-5 && k<20']) % ��n��������������ֹ����
        k=k+1;
        K_Global=zeros(N_Node*3,N_Node*3); % �γɿյ��ܸվ���
        for i=1:N_Element    % �����ܸ�
            N1=Element(i,1); % ��ȡ��Ԫi�ڵ���
            N2=Element(i,2); % ��ȡ��Ԫj�ڵ���
            C(i,:)=Update_Node(N2,:)-Update_Node(N1,:); % �������ɵĵ�Ԫ��������
            L(i)=norm(C(i,:)); % �������ɵĵ�Ԫ���ȱ仯
            cs=C(i,:)/L(i);    % ��Ԫת�ǵ�cos��sinֵ
            lamda=lamdabianhuan(cs); % ����������ֲ�����ת������
            K_Element_jb=beam2d_stiffness(i,L(i),Element_F(i,:),Element,Section); % �õ���Ԫ�ھֲ������µĸնȾ���
            K_Element_zt=lamda'*K_Element_jb*lamda; % �õ���Ԫ�����������µĸնȾ���
            K_Global=K_Global+Assemble(K_Element_zt,Element(i,:),N_Node); % ��װ���յõ��ܸ�
        end
        eval(['[XZK,XZP]=bianjiexiuzheng(Boundary,Val',num2str(n),'_',num2str(k-1),',K_Global,N_Node);']) % �ں��زв������µı߽�����
        dis2=XZK\XZP; % �����زв�µ�λ��
        [XZK,XZP]=bianjiexiuzheng(Boundary,Force_Transfer,K_Global,N_Node); % ���б߽�����
        dis3=XZK\XZP; % �������ص��µ�λ��
        lamda=-dis1'*dis2/(dis1'*dis3+1); % ������ϵ��
        Disp_Transfer=lamda*dis3+dis2; % �������λ��
        r_k=[Disp_Transfer',lamda]'; % ������¼
        R_m=R_m+r_k; % ��������
        zlamda=zlamda+lamda; % �ܺ���ϵ������
        All_F=zlamda*Force_Transfer; % �ܺ��ظ���
        for i=1:N_Node
            All_Disp(i,:)=All_Disp(i,:)+Disp_Transfer(3*i-2:3*i)'; % ��λ�Ƹ���
        end
        for i=1:N_Element;
            Disp_Element_zt=zeros(6,1); % ���������µ�Ԫλ�ƴ洢����
            N1=Element(i,1); % ��ȡ��Ԫi�ڵ���
            N2=Element(i,2); % ��ȡ��Ԫj�ڵ���
            for j=1:3 % ��ȡ��Ԫi�ڵ�λ��
                Disp_Element_zt(j,1)=Disp_Transfer(3*N1+j-3);
            end
            for j=4:6 % ��ȡ��Ԫi�ڵ�λ��
                Disp_Element_zt(j,1)=Disp_Transfer(3*N2+j-6);
            end
            C(i,:)=Update_Node(N2,:)-Update_Node(N1,:); % �������ɵĵ�Ԫ��������
            L(i)=norm(C(i,:)); % �������ɵĵ�Ԫ���ȱ仯
            cs=C(i,:)/L(i); % ��Ԫת�ǵ�cos��sinֵ
            lamda=lamdabianhuan(cs); % ����������ֲ�����ת������
            Disp_Element_jb=lamda*Disp_Element_zt; % �õ���Ԫ�ھֲ������µı���
            Disp_Element_jb_tanxing=tanxingbianxing(Disp_Element_jb,L(i)); % �õ���Ԫ�ھֲ������µĵ��Ա���
            K_Element_jb=beam2d_stiffness(i,L(i),Element_F(i,:),Element,Section); % �õ���Ԫ�ھֲ������µĸնȾ���
            Element_F_jb=K_Element_jb*Disp_Element_jb_tanxing; % �õ���Ԫ�ھֲ������µ�����
            Element_F(i,:)=Element_F(i,:)+Element_F_jb(:,1)';  % ���µ�Ԫ����
        end
        Update_Node=Node+All_Disp(:,1:2); % ����λ����Ϣ���µĽڵ�����,x������y����
        Internal_F=zeros(N_Node*3,1); % �յĽڵ������洢����
        for i=1:N_Element
            F_Global=zeros(N_Node*3,1); % �������������½ڵ��������
            N1=Element(i,1); % ��ȡ��Ԫi�ڵ���
            N2=Element(i,2); % ��ȡ��Ԫj�ڵ���
            C(i,:)=Update_Node(N2,:)-Update_Node(N1,:); % �������ɵĵ�Ԫ��������
            L(i)=norm(C(i,:)); % �������ɵĵ�Ԫ���ȱ仯
            cs=C(i,:)/L(i); % ��Ԫת�ǵ�cos��sinֵ
            lamda=lamdabianhuan(cs); % ����������ֲ�����ת������
            Element_All_F_jb=zeros(6,1); % ���ɾֲ��ڵ���������洢����
            Element_All_F_jb(:,1)=Element_F(i,:); % ���ɾֲ��ڵ��������
            Element_All_F_zt=lamda'*Element_All_F_jb; % �õ���������ϵ�½ڵ����
            F_Global(3*N1-2:3*N1,1)=Element_All_F_zt(1:3,1); % �õ����������½ڵ����
            F_Global(3*N2-2:3*N2,1)=Element_All_F_zt(4:6,1); % �õ����������½ڵ����
            Internal_F=Internal_F+F_Global; % �ڵ���������
        end
        eval(['Val',num2str(n),'_',num2str(k),'=All_F-Internal_F;']) % �õ���n����������k�ε����ĺ��زв�
        for i=1:N_Node % �Ժ��زв���б߽�����
            for j=2:4
                if Boundary(i,j)==0
                   eval(['Val',num2str(n),'_',num2str(k),'(3*i+j-4)=0;'])
                end
            end
        end
   end
   Y=[Y,0];
   X=[X,0];
   Y(1,n+1)=-All_F(17,1);   % ��¼���к���
   X(1,n+1)=-All_Disp(6,2); % ��¼��������λ��
end
% ��ͼ
plot(X,Y);
grid on;