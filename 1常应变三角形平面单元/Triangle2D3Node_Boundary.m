function [z,p]=Triangle2D3Node_Boundary(KK,P,DOF,disp)
% �ú����ó˴������޸��ܸնȾ���ʵ�ֱ߽�����
% �����ܸնȾ���KK����������P��ҪԼ�������ɶ�DOF��Լ��λ��ֵdisp
% ����ܸնȾ���z����������p

BigNumber=10^20;
P(DOF,1)=BigNumber*KK(DOF,DOF)*disp;
p=P;
KK(DOF,DOF)=KK(DOF,DOF)*BigNumber;
z=KK;
end