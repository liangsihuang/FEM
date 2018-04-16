function stress=Triangle2D3Node_Stress(ie,ID)
% �ú������㵥Ԫ��Ӧ��
% �����Ԫ��Ӧ��stress(3��1),������Ϊ��Ӧ����Ԫ����Ԫ��Ӧ������ΪSx,Sy,Sxy
global gNode gElement gMaterial gDelta
% ��ȡ�ڵ����꣺�����ڵ�i��j��m������xi,yi,xj,yj,xm,ym
xi=gNode(gElement(ie,1),1);
yi=gNode(gElement(ie,1),2);
xj=gNode(gElement(ie,2),1);
yj=gNode(gElement(ie,2),2);
xm=gNode(gElement(ie,3),1);
ym=gNode(gElement(ie,3),2);
% ��ȡ���t
t=gElement(ie,5);
% ��ȡ����ģ��E�����ɱ�NU
E=gMaterial(gElement(ie,4),1);
NU=gMaterial(gElement(ie,4),2);
% ��������
ai=xj*ym-xm*yj;aj=xm*yi-xi*ym;am=xi*yj-xj*yi;
bi=yj-ym;bj=ym-yi;bm=yi-yj;
ci=-xj+xm;cj=-xm+xi;cm=-xi+xj;
%���������
A=1/2*(ai+aj+am);
% ���ξ���
B=1/(2*A)*[bi 0 bj 0 bm 0;
    0 ci 0 cj 0 cm;
    ci bi cj bj cm bm];

if ID==1
    D=E/(1-NU*NU)*[1 NU 0;NU 1 0;0 0 (1-NU)/2];
elseif ID==2
    D=E/(1-2*NU)/(1+NU)*[1-NU NU 0;NU 1-NU 0;0 0 (1-2*NU)/2];
end
% ��ȡ��Ԫ��λ������u��6��1��
u=zeros(6,1);
u(1,1)=gDelta((gElement(ie,1)-1)*2+1);
u(2,1)=gDelta((gElement(ie,1)-1)*2+2);
u(3,1)=gDelta((gElement(ie,2)-1)*2+1);
u(4,1)=gDelta((gElement(ie,2)-1)*2+2);
u(5,1)=gDelta((gElement(ie,3)-1)*2+1);
u(6,1)=gDelta((gElement(ie,3)-1)*2+2);

stress=D*B*u;
end