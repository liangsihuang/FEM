function BernoulliBeam2D_initialModel
%  定义平面杆系的有限元模型，单位N，m
%        gNode ------- 节点定义
%        gElement ---- 单元定义
%        gMaterial --- 材料定义，包括弹性模量，梁的截面积和梁的抗弯惯性矩
%        gBC1 -------- 约束条件
%        gNF --------- 集中力
%        gDF --------- 分布力

    global gNode gElement gMaterial gBC1 gNF gDF
    
[~,SheetNames]=xlsfinfo('model.xlsx');
nSheets=length(SheetNames);
for i = 1:nSheets
  name = SheetNames{i}; 
  data = readtable('model.xlsx','Sheet',name) ; 
  S(i).name = name;
  S(i).data = data;
end

    % 节点坐标 x  y
  gNode(:,1) = S(1).data.X;
  gNode(:,2) = S(1).data.Y;
    % 单元定义  节点1  节点2  材料号
  gElement(:,1) = S(2).data.Node1;
  gElement(:,2) = S(2).data.Node2;
  gElement(:,3) = S(2).data.Material;  
    % 材料性质  弹性模量  抗弯惯性矩  截面积
    gMaterial = [2.06e11 1.0e-4 1.0e-2];

    % 第一类约束条件 节点号   自由度号    约束值
    gBC1 = [ 1 1 0.0;
             1 2 0.0;
             1 3 0.0;];
         
E=2.06*10^11; %Pa
I=1*10^(-4); %m^3
L=1; %m
M=E*I*(2*pi/L); %N*m
% 集中力 节点号   自由度号   集中力值
gNF = [11 3 -0.01*M];

return