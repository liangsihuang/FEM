clc;clear;
BuildModel();
expectU=-2;
record=SolveModel(expectU);
global Node Element Material BC1 NF TotalU InterF_elem InterF ExterF dU r_total lambda_total
x=record(:,1)*(-1);
y=record(:,2);
plot(x,y);
% 'o'