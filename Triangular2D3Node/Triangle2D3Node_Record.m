function Triangle2D3Node_Record
%  ��¼������
global gNode gElement gMaterial gBC1 gNF gDF gK gDelta gStress
    
fid=fopen('node_disp.txt','w');
[node_number,dummy] = size( gNode ) ;
for i=1:node_number
%     ΢���Դ�txt�ı���������ʹ��\r\n��һ��ʹ��\n
    fprintf(fid,'%i %f %f\r\n',i,gDelta((i-1)*2+1),gDelta((i-1)*2+2)); 
end
fclose(fid);

fid=fopen('element_stress.txt','w');
[element_number,dummy] = size( gElement ) ;
for i=1:element_number
%     ΢���Դ�txt�ı���������ʹ��\r\n��һ��ʹ��\n
    fprintf(fid,'%i %f %f %f\r\n',i,gStress(i,1),gStress(i,2),gStress(i,3)); 
end
fclose(fid);
return