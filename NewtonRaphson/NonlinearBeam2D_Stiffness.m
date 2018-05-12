function k=NonlinearBeam2D_Stiffness(ie,icoord)
% 输入:
% ie:单元号
% icoord:坐标系参数，1代表整体坐标系，2代表局部坐标系
% 输出：
% k：根据icoord的值，相应坐标系下的单元刚度矩阵
global gNode gElement gMaterial
k = zeros( 6, 6 ) ;
E = gMaterial( gElement(ie, 3), 1 ) ;
I = gMaterial( gElement(ie, 3), 2 ) ;
A = gMaterial( gElement(ie, 3), 3 ) ;
xi = gNode( gElement( ie, 1 ), 1 ) ;
yi = gNode( gElement( ie, 1 ), 2 ) ;
xj = gNode( gElement( ie, 2 ), 1 ) ;
yj = gNode( gElement( ie, 2 ), 2 ) ;
L = ( (xj-xi)^2 + (yj-yi)^2 )^(1/2) ;
k_elastic = [  E*A/L           0          0 -E*A/L           0          0
                   0  12*E*I/L^3  6*E*I/L^2      0 -12*E*I/L^3  6*E*I/L^2
                   0   6*E*I/L^2    4*E*I/L      0  -6*E*I/L^2    2*E*I/L
              -E*A/L           0          0  E*A/L           0          0
                   0 -12*E*I/L^3 -6*E*I/L^2      0  12*E*I/L^3 -6*E*I/L^2
                   0   6*E*I/L^2    2*E*I/L      0  -6*E*I/L^2    4*E*I/L];
enf=NonlinearBeam2D_ElementNodeForce(ie);
F=enf(1,4); %为什么不用左端点的轴力值？跟正负号有关
% 这里有一个假定，没有轴向分布力
M1=enf(1,3);
M2=enf(1,6);
k_geom = [  F/L                        0                 -M1/L   -F/L                        0                -M2/L
              0             12*F*I/A/L^3      6*F*I/A/L^2+F/10      0   -(12*F*I/A/L^3+6*F/5/L)    6*F*I/A/L^2+F/10
          -M1/L         6*F*I/A/L^2+F/10    4*F*I/A/L+2*F*L/15   M1/L       -(6*F*I/A/L^2+F/10)    2*F*I/A/L-F*L/30
           -F/L                        0                  M1/L    F/L                        0                 M2/L
              0   -(12*F*I/A/L^3+6*F/5/L)   -(6*F*I/A/L^2+F/10)     0      12*F*I/A/L^3+6*F/5/L  -(6*F*I/A/L^2+F/10)
          -M2/L         6*F*I/A/L^2+F/10      2*F*I/A/L-F*L/30   M2/L        -(6*F*I/A/L^2+F/10)   4*F*I/A/L+2*F*L/15];
k=k_elastic+k_geom;
if icoord == 1
    T = NonlinearBeam2D_TransformMatrix( ie );
    k = T*k*transpose(T);
end
end