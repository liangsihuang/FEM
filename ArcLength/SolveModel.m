function record=SolveModel()
contU=0;
increNum=0;
record=zeros(1000,2);
while(contU<-2) %���㵽����λ��Ϊ2m
    increNum=increNum+1;
    [lambda ,r_total]=ArcLengthMethod(increNum, lambda, r_total);
    if (increNum>size(record,1))
        temp=zeros(1000,2);
        record=[record,temp]; %ֻ�ں�������Ż�ı�record��ά��,���ۿ��Խ���
    end
    contU=MidSpanU(6,2);
    record(increNum,1)=contU;
    record(increNum,2)=lambda;
end
record=record(increNum,:);
end

function u=MidSpanU(n,d)
% n�����ƵĽڵ��
% d�����ɶ�
% u���ڵ�λ��
global TotalU
u=TotalU(3*(n-1)+d);
end

