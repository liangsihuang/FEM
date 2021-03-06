function record=SolveModel(expectU)
global Node r_total lambda_total dr_total
contU=0;
increNum=0;
record=zeros(1000,2);
lambda_total=0;
r_total=zeros(3*size(Node,1)+1,1);
dr_total=zeros(3*size(Node,1)+1,1); %m-1到m的向量（delta弧长）
while(contU>expectU) 
    increNum=increNum+1;
    ArcLengthMethod(increNum);
    if (increNum>size(record,1))
        temp=zeros(1000,2);
        record=[record;temp]; %只在很少情况才会改变record的维度,代价可以接受
    end
    contU=ControlPointU();
    record(increNum,1)=contU;
    record(increNum,2)=lambda_total;
end
record=record(1:increNum,:);
end

function u=ControlPointU() % 控制点为加载点
% n：控制的节点号
% d：自由度
% u：节点位移
global TotalU NF
n = NF( 1, 1 );
d = NF( 1, 2 );
u=TotalU(3*(n-1)+d);
end

