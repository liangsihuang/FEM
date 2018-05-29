function UpdateNode()
global Node dU
[num,~]=size(Node);
for i=1:num
   Node(i,1)=Node(i,1)+dU(i*3-2);
   Node(i,2)=Node(i,2)+dU(i*3-1);
end
end