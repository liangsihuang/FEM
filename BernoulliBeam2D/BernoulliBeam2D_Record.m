function BernoulliBeam2D_Record
%  ��¼������
global gNode gElement gMaterial gBC1 gNF gDF gK gDelta gStress
    
fid=fopen('node_disp.txt','w');
[node_number,dummy] = size( gNode ) ;
for i=1:node_number
%     ΢���Դ�txt�ı���������ʹ��\r\n��һ��ʹ��\n
    fprintf(fid,'%i %f %f %f\r\n',i,gDelta((i-1)*3+1)+0,gDelta((i-1)*3+2)+0,gDelta((i-1)*3+1)+0); 
%     ����0�ᱨ������ʹ�� fprintf��û��Ϊϡ�����붨�庯����
end
fclose(fid);

fid=fopen('element_node_force.txt','w');
[element_number, dummy] = size( gElement ) ;
for ie = 1:element_number
    enf = BernoulliBeam2D_ElementNodeForce( ie );
%     ΢���Դ�txt�ı���������ʹ��\r\n��һ��ʹ��\n
    fprintf(fid,'%i %i %f %f %f\r\n',ie,gElement(ie,1),enf(1),enf(2),enf(3)); 
    fprintf(fid,'%i %i %f %f %f\r\n',ie,gElement(ie,2),enf(4),enf(5),enf(6));
%     �����У���һ��Ϊ��Ԫ�ĵ�һ���ڵ㴦�����������������
end
fclose(fid);

return