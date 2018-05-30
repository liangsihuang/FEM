function TestAssembleForce()
global InterF f_onlyNF InterF_elem Element
[m,n]=size(f_onlyNF);
InterF=zeros(m,n);
CreateInterF_elem();
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

function CreateInterF_elem()
global Element InterF_elem
[num,~]=size(Element);
for ie=1:num
%     enf=ElementInternalForce(ie);
    enf = BernoulliBeam2D_ElementNodeForce( ie );
    InterF_elem(:,ie)=enf;
end

end