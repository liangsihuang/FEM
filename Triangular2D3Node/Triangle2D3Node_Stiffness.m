function k=Triangle2D3Node_Stiffness(ie,ID)
% 该函数计算单元的刚度矩阵
% 输入单元号
% 输入平面问题性质指标参数ID（1为平面应力，2为平面应变）
% 输出单元刚度矩阵k（6×6）
global gNode gElement gMaterial
% 提取节点坐标：三个节点i、j、m的坐标xi,yi,xj,yj,xm,ym
xi=gNode(gElement(ie,1),1);
yi=gNode(gElement(ie,1),2);
xj=gNode(gElement(ie,2),1);
yj=gNode(gElement(ie,2),2);
xm=gNode(gElement(ie,3),1);
ym=gNode(gElement(ie,3),2);
% 提取厚度t
t=gElement(ie,5);
% 提取弹性模量E，泊松比NU
E=gMaterial(gElement(ie,4),1);
NU=gMaterial(gElement(ie,4),2);

% 参数计算
ai=xj*ym-xm*yj;aj=xm*yi-xi*ym;am=xi*yj-xj*yi;
bi=yj-ym;bj=ym-yi;bm=yi-yj;
ci=-xj+xm;cj=-xm+xi;cm=-xi+xj;
%三角形面积，取绝对值就不管顺时针逆时针了？
A=abs(1/2*(ai+aj+am));
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