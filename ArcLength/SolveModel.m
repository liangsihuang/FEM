function record=SolveModel()
global Node
contU=0;
increNum=0;
record=zeros(1000,2);
lambda_total=0;
r_total=zeros(3*size(Node,1)+1,1);
while(contU>-2) %���㵽����λ��Ϊ2m
    increNum=increNum+1;
    [lambda_total ,r_total]=ArcLengthMethod(increNum, lambda_total, r_total);
    if (increNum>size(record,1))
        temp=zeros(1000,2);
        record=[record;temp]; %ֻ�ں�������Ż�ı�record��ά��,���ۿ��Խ���
    end
    contU=ControlPointU(6,2);
    record(increNum,1)=contU;
    record(increNum,2)=lambda_total;
end
record=record(increNum,:);
end

function u=ControlPointU(n,d)
% n�����ƵĽڵ��
% d�����ɶ�
% u���ڵ�λ��
global TotalU
u=TotalU(3*(n-1)+d);
end

