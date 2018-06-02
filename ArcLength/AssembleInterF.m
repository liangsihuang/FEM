function AssembleInterF()
global InterF_elem InterF Element
% ±ÿ–Îœ»«Â¡„
[m,n]=size(InterF);
InterF=zeros(m,n);
[~,num]=size(InterF_elem);
for ie=1:num
    n1=Element(ie,1);
    n2=Element(ie,2);
    T=Beam2D_TransformMatrix(ie);
    temp=T'*InterF_elem(:,ie);
    InterF(n1*3-2:n1*3)=InterF(n1*3-2:n1*3)+temp(1:3);
    InterF(n2*3-2:n2*3)=InterF(n2*3-2:n2*3)+temp(4:6);
end

end
