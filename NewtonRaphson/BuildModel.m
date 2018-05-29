function BuildModel()
global Node Element Material BC1 NF Disp InterF_elem InterF ExterF
% InterF_elem: 单元内力矩阵，行数=6，列数=单元数
[~,SheetNames]=xlsfinfo('model.xlsx');
nSheets=length(SheetNames);
for i = 1:nSheets
  name = SheetNames{i};
  data = readtable('model.xlsx','Sheet',name) ; 
  S(i).name = name;
  S(i).data = data;
end
% 导入的数据
Node(:,1) = S(1).data.X;
Node(:,2) = S(1).data.Y;

Element(:,1) = S(2).data.Node1;
Element(:,2) = S(2).data.Node2;
Element(:,3) = S(2).data.Material;

Material(:,1) = S(3).data.E;
Material(:,2) = S(3).data.I;
Material(:,3) = S(3).data.A;

BC1(:,1) = S(4).data.Node;
BC1(:,2) = S(4).data.DOF;
BC1(:,3) = S(4).data.Constrain;

NF(:,1) = S(5).data.Node;
NF(:,2) = S(5).data.DOF;
NF(:,3) = S(5).data.Force;

% 初始化为0的数据
[m,~]=size(Node);
Disp=zeros(m*3,1);
InterF=zeros(m*3,1);
ExterF=zeros(m*3,1);
% K=zeros(m*3,m*3);

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