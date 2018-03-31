function stress=Triangle2D3Node_Stress(E,NU,xi,yi,xj,yj,xm,ym,u,ID)
% �ú������㵥Ԫ��Ӧ��
% ���뵯��ģ��E,���ɱ�NU�����t
% ���������ڵ�i,j,m������xi,yi,xj,yj,xm,ym
% ����ƽ����������ָ�����ID��1Ϊƽ��Ӧ����2Ϊƽ��Ӧ�䣩
% ���뵥Ԫ��λ������u��6��1��
% �����Ԫ��Ӧ��stress(3��1),������Ϊ��Ӧ����Ԫ����Ԫ��Ӧ������ΪSx,Sy,Sxy

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

stress=D*B*u;
end