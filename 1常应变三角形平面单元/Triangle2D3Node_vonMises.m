function von=Triangle2D3Node_vonMises(stress)
% 将应力向量（3*1）转化为vonMises应力值
x=stress(1,1);
y=stress(2,1);
t=stress(3,1);
sigma_1=1/2*(x+y)+sqrt(((x-y)/2)^2+t^2);
sigma_2=1/2*(x+y)-sqrt(((x-y)/2)^2+t^2);
von=sqrt(sigma_1^2+sigma_2^2-sigma_1*sigma_2);
end