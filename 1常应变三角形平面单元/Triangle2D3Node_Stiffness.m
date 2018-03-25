function k=Triangle2D3Node_Stiffness(E,NU,t,xi,yi,xj,yj,xm,ym,ID)
% �ú������㵥Ԫ�ĸնȾ���
% ���뵯��ģ��E�����ɱ�NU�����t
% ���������ڵ�i��j��m������xi,yi,xj,yj,xm,ym
% ����ƽ����������ָ�����ID��1Ϊƽ��Ӧ����2Ϊƽ��Ӧ�䣩
% �����Ԫ�նȾ���k��6��6��

%���������
A=1/2*(xi*(yj-ym)+xj*(ym-yi)+xm*(yi-yj));
% ��������
bi=yj-ym;bj=ym-yi;bm=yi-yj;
ci=-xj+xm;cj=-xm+xi;cm=-xi+xj;
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