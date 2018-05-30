function UpdateInterF()
% ���� InterF (��������)
global InterF InterF_elem Element
[num,~]=size(Element);
for ie=1:num
    n1=Element(ie,1);
    n2=Element(ie,2);
    d_inf=InterF_elem(:,ie);
    T=Beam2D_TransformMatrix(ie); %������C2״̬�µ�ת������Ҫ��updateNode
    g_d_inf=T'*d_inf;
    InterF(n1*3-2:n1*3)=InterF(n1*3-2:n1*3)+g_d_inf(1:3);
    InterF(n2*3-2:n2*3)=InterF(n2*3-2:n2*3)+g_d_inf(4:6);
end

end



