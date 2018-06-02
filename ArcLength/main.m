clc;clear;
BuildModel();
record=SolveModel();
global Node Element Material BC1 NF TotalU InterF_elem InterF ExterF dU
Plot(record);