function BernoulliBeam2D_BuildModel()
% ����ƽ���ϵ������Ԫģ�ͣ���λN��m
% Node ------- �ڵ㶨��
% Element ---- ��Ԫ����
% Material --- ���϶��壬��������ģ�������Ľ���������Ŀ�����Ծ�
% BC1 -------- Լ������
% NF --------- ������
% DF --------- �ֲ���

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