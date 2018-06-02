function BuildModel()
global Node Element Material BC1 NF TotalU InterF_elem InterF ExterF
% InterF_elem: 单元内力矩阵，行数=6，列数=单元数
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
end

% 初始化为0的数据
[m,~]=size(Node);
InterF=zeros(m*3,1);
ExterF=zeros(m*3,1);
TotalU=zeros(m*3,1);

[n,~]=size(Element);
InterF_elem=zeros(6,n);


% 把集中力直接集成到整体节点力向量中
[num, ~] = size( NF ) ;
for i=1:1:num
    n = NF( i, 1 ) ;
    d = NF( i, 2 ) ;
    ExterF( (n-1)*3 + d ) = NF( i, 3 ) ;
end


end