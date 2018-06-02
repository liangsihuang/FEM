function record=SolveModel()
contU=0;
increNum=0;
record=zeros(1000,2);
while(contU<-2) %计算到跨中位移为2m
    increNum=increNum+1;
    [lambda ,r_total]=ArcLengthMethod(increNum, lambda, r_total);
    if (increNum>size(record,1))
        temp=zeros(1000,2);
        record=[record,temp]; %只在很少情况才会改变record的维度,代价可以接受
    end
    contU=MidSpanU(6,2);
    record(increNum,1)=contU;
    record(increNum,2)=lambda;
end
record=record(increNum,:);
end

function u=MidSpanU(n,d)
% n：控制的节点号
% d：自由度
% u：节点位移
global TotalU
u=TotalU(3*(n-1)+d);
end

