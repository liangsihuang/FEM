function K_G=Assemble(K_Element,Element,N_Node)
% 输入：
% K_Element:单元刚度矩阵
% Element:
% N_Node:
% 输出：
% K_G: global stiffness matrix
K_G=zeros(N_Node*3,N_Node*3); %???
for i=1:2
    n1=Element(1,i);
    for j=1:2
        n2=Element(1,j);
        K_G(3*n1-2:3*n1,3*n2-2:3*n2)=K_Element(3*i-2:3*i,3*j-2:3*j);
    end
end
end


