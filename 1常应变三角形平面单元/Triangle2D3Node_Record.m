function Triangle2D3Node_Record
%  ��¼������
global gNode gElement gMaterial gBC1 gNF gDF gK gDelta
    
fid=fopen('node0.txt','w');
[node_number,dummy] = size( gNode ) ;
for i=1:node_number
%     ΢���Դ�txt�ı���������ʹ��\r\n��һ��ʹ��\n
    fprintf(fid,'%i %f %f\r\n',i,gDelta((i-1)*2+1),gDelta((i-1)*2+2)); 
end
fclose(fid);
return