function BernoulliBeam2D_initialModel
%  ����ƽ���ϵ������Ԫģ�ͣ���λN��m
%        gNode ------- �ڵ㶨��
%        gElement ---- ��Ԫ����
%        gMaterial --- ���϶��壬��������ģ�������Ľ���������Ŀ�����Ծ�
%        gBC1 -------- Լ������
%        gNF --------- ������
%        gDF --------- �ֲ���

    global gNode gElement gMaterial gBC1 gNF gDF
    
[~,SheetNames]=xlsfinfo('model.xlsx');
nSheets=length(SheetNames);
for i = 1:nSheets
  name = SheetNames{i}; 
  data = readtable('model.xlsx','Sheet',name) ; 
  S(i).name = name;
  S(i).data = data;
end

    % �ڵ����� x  y
  gNode(:,1) = S(1).data.X;
  gNode(:,2) = S(1).data.Y;
    % ��Ԫ����  �ڵ�1  �ڵ�2  ���Ϻ�
  gElement(:,1) = S(2).data.Node1;
  gElement(:,2) = S(2).data.Node2;
  gElement(:,3) = S(2).data.Material;  
    % ��������  ����ģ��  ������Ծ�  �����
    gMaterial = [2.06e11 1.0e-4 1.0e-2];

    % ��һ��Լ������ �ڵ��   ���ɶȺ�    Լ��ֵ
    gBC1 = [ 1 1 0.0;
             1 2 0.0;
             1 3 0.0;];
         
E=2.06*10^11; %Pa
I=1*10^(-4); %m^3
L=1; %m
M=E*I*(2*pi/L); %N*m
% ������ �ڵ��   ���ɶȺ�   ������ֵ
gNF = [11 3 -0.01*M];

return