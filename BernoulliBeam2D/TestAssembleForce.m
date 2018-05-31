function TestAssembleForce()
global InterF f_onlyNF InterF_elem Element R
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

R=InterF-f_onlyNF;
R=ModifyResidual(R);
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

function R=ModifyResidual(R)
global BC1
[num,~] = size( BC1 ) ;
for i=1:1:num
    n = BC1(i, 1);
    d = BC1(i, 2);
    m = (n-1)*3 + d;
    R(m) = BC1(i, 3);
end
end