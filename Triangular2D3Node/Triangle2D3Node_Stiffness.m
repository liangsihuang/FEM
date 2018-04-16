function k=Triangle2D3Node_Stiffness(ie,ID)
% �ú������㵥Ԫ�ĸնȾ���
% ���뵥Ԫ��
% ����ƽ����������ָ�����ID��1Ϊƽ��Ӧ����2Ϊƽ��Ӧ�䣩
% �����Ԫ�նȾ���k��6��6��
global gNode gElement gMaterial
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
%�����������ȡ����ֵ�Ͳ���˳ʱ����ʱ���ˣ�
A=abs(1/2*(ai+aj+am));
% ���ξ���
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