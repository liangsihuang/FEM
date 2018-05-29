function BuildModel()
global Node Element Material BC1 NF Disp
[~,SheetNames]=xlsfinfo('model.xlsx');
nSheets=length(SheetNames);
for i = 1:nSheets
  name = SheetNames{i}; 
  data = readtable('model.xlsx','Sheet',name) ; 
  S(i).name = name;
  S(i).data = data;
end

Node(:,1) = S(1).data.X;
Node(:,2) = S(1).data.Y;

Element(:,1) = S(2).data.Node1;
Element(:,2) = S(2).data.Node2;
Element(:,3) = S(2).data.Material;

Material(:,1) = S(3).data.E;
Material(:,2) = S(3).data.I;
Material(:,3) = S(3).data.A;

% 第一类约束条件 节点号   自由度号    约束值
BC1(:,1) = S(4).data.Node;
BC1(:,2) = S(4).data.DOF;
BC1(:,3) = S(4).data.Constrain;

NF(:,1) = S(5).data.Node;
NF(:,2) = S(5).data.DOF;
NF(:,3) = S(5).data.Force;

[m,n]=sizes(Node);
Disp=zeros(m,n);
end