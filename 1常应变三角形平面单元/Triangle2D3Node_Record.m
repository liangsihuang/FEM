function Triangle2D3Node_Record
%  记录计算结果
global gNode gElement gMaterial gBC1 gNF gDF gK gDelta
    
fid=fopen('node0.txt','w');
[node_number,dummy] = size( gNode ) ;
for i=1:node_number
%     微软自带txt文本打开器换行使用\r\n，一般使用\n
    fprintf(fid,'%i %f %f\r\n',i,gDelta((i-1)*2+1),gDelta((i-1)*2+2)); 
end
fclose(fid);
return