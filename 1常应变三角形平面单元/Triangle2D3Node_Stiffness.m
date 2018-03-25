function k=Triangle2D3Node_Stiffness(E,NU,t,xi,yi,xj,yj,xm,ym,ID)
% 该函数计算单元的刚度矩阵
% 输入弹性模量E，泊松比NU，厚度t
% 输入三个节点i、j、m的坐标xi,yi,xj,yj,xm,ym
% 输入平面问题性质指标参数ID（1为平面应力，2为平面应变）
% 输出单元刚度矩阵k（6×6）

%三角形面积
A=1/2*(xi*(yj-ym)+xj*(ym-yi)+xm*(yi-yj));
% 参数计算
bi=yj-ym;bj=ym-yi;bm=yi-yj;
ci=-xj+xm;cj=-xm+xi;cm=-xi+xj;
% 几何矩阵
B=1/(2*A)*[bi 0 bj 0 bm 0;
    0 ci 0 cj 0 cm;
    ci bi cj bj cm bm];
if ID==1
    D=E/(1-NU*NU)*[1 NU 0;NU 1 0;0 0 (1-NU)/2];
elseif ID==2
    D=E/(1-2*NU)/(1+NU)*[1-NU NU 0;NU 1-NU 0;0 0 (1-2*NU)/2];
end

k=B'*D*B*t*A;
end