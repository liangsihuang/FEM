function BernoulliBeam2D_BuildModel()
% 定义平面杆系的有限元模型，单位N，m
% Node ------- 节点定义
% Element ---- 单元定义
% Material --- 材料定义，包括弹性模量，梁的截面积和梁的抗弯惯性矩
% BC1 -------- 约束条件
% NF --------- 集中力
% DF --------- 分布力

global Node Element Material BC1 NF DF

[~,sheets]=xlsfinfo('model.xlsx');
nSheets=length(sheets);
for i = 1:nSheets
    sheet=cell2mat(sheets(i));
    if (strcmp(sheet,'node'))
        Node=xlsread('model.xlsx',sheet);
    end
    if (strcmp(sheet,'element'))
        Element=xlsread('model.xlsx',sheet);
    end
    if (strcmp(sheet,'material'))
        Material=xlsread('model.xlsx',sheet);
    end
    if (strcmp(sheet,'boundary'))
        BC1=xlsread('model.xlsx',sheet);
    end
    if (strcmp(sheet,'nodalforce'))
        NF=xlsread('model.xlsx',sheet);
    end
    if (strcmp(sheet,'distributedforce'))
        DF=xlsread('model.xlsx',sheet);
    end
end
return