function k=NonlinearBeam2D_Stiffness(ie,icoord, inf)
% 输入:
% ie:单元号
% icoord:坐标系参数，1代表整体坐标系，2代表局部坐标系
% inf: initial(internal) force, 单元初应力(内力)
% 输出：
% k：根据icoord的值，相应坐标系下的单元刚度矩阵
global Node Element Material
k = zeros( 6, 6 ) ;
E = Material( Element(ie, 3), 1 ) ;
I = Material( Element(ie, 3), 2 ) ;
A = Material( Element(ie, 3), 3 ) ;
xi = Node( Element( ie, 1 ), 1 ) ;
yi = Node( Element( ie, 1 ), 2 ) ;
xj = Node( Element( ie, 2 ), 1 ) ;
yj = Node( Element( ie, 2 ), 2 ) ;
L = ( (xj-xi)^2 + (yj-yi)^2 )^(1/2) ;
ke = [  E*A/L           0          0 -E*A/L           0          0
                   0  12*E*I/L^3  6*E*I/L^2      0 -12*E*I/L^3  6*E*I/L^2
                   0   6*E*I/L^2    4*E*I/L      0  -6*E*I/L^2    2*E*I/L
              -E*A/L           0          0  E*A/L           0          0
                   0 -12*E*I/L^3 -6*E*I/L^2      0  12*E*I/L^3 -6*E*I/L^2
                   0   6*E*I/L^2    2*E*I/L      0  -6*E*I/L^2    4*E*I/L];
F=inf(4,1); %为什么不用左端点的轴力值？跟正负号有关。这里有一个假定，没有轴向分布力
M1=inf(3,1);
M2=inf(6,1);
kg =     [  F/L                        0                 -M1/L   -F/L                        0                -M2/L
              0             12*F*I/A/L^3      6*F*I/A/L^2+F/10      0   -(12*F*I/A/L^3+6*F/5/L)    6*F*I/A/L^2+F/10
          -M1/L         6*F*I/A/L^2+F/10    4*F*I/A/L+2*F*L/15   M1/L       -(6*F*I/A/L^2+F/10)    2*F*I/A/L-F*L/30
           -F/L                        0                  M1/L    F/L                        0                 M2/L
              0   -(12*F*I/A/L^3+6*F/5/L)   -(6*F*I/A/L^2+F/10)     0      12*F*I/A/L^3+6*F/5/L  -(6*F*I/A/L^2+F/10)
          -M2/L         6*F*I/A/L^2+F/10      2*F*I/A/L-F*L/30   M2/L        -(6*F*I/A/L^2+F/10)   4*F*I/A/L+2*F*L/15];
k=ke+kg;
if icoord == 1
    T = Beam2D_TransformMatrix( ie );
    k = T*k*transpose(T);
end
end