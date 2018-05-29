function AssembleInterF()
global InterF_elem InterF Element
[~,num]=sizes(InterF_elem);
for ie=1:num
    n1=Element(ie,1);
    n2=Element(ie,2);
    InterF(n1*3-2:n1*3)=InterF(n1*3-2:n1*3)+InterF_elem(1:3,ie);
    InterF(n2*3-2:n2*3)=InterF(n2*3-2:n2*3)+InterF_elem(4:6,ie);
end

end
