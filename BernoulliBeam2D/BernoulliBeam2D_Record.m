function BernoulliBeam2D_Record
%  ��¼������
global Node Element dU
    
fid=fopen('node_disp.txt','w');
[node_number,~] = size( Node ) ;
for i=1:node_number
%     ΢���Դ�txt�ı���������ʹ��\r\n��һ��ʹ��\n
    fprintf(fid,'%i %f %f %f\r\n',i,dU((i-1)*3+1)+0,dU((i-1)*3+2)+0,dU((i-1)*3+1)+0); 
%     ����0�ᱨ������ʹ�� fprintf��û��Ϊϡ�����붨�庯����
end
fclose(fid);

fid=fopen('element_node_force.txt','w');
[element_number, ~] = size( Element ) ;
for ie = 1:element_number
    enf = BernoulliBeam2D_ElementNodeForce( ie );
%     ΢���Դ�txt�ı���������ʹ��\r\n��һ��ʹ��\n
    fprintf(fid,'%i %i %f %f %f\r\n',ie,Element(ie,1),enf(1),enf(2),enf(3)); 
    fprintf(fid,'%i %i %f %f %f\r\n',ie,Element(ie,2),enf(4),enf(5),enf(6));
%     �����У���һ��Ϊ��Ԫ�ĵ�һ���ڵ㴦�����������������
    
end
fclose(fid);

return