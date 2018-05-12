clc;clear;
E=2.06*10^11; %Pa
I=1*10^(-4); %m^3
L=1; %m
M=E*I*(2*pi/L); %N*m

[~,SheetNames]=xlsfinfo('model.xlsx');
nSheets=length(SheetNames);
for i = 1:nSheets
  name = SheetNames{i}; 
  data = readtable('model.xlsx','Sheet',name) ; 
  S(i).name = name;
  S(i).data = data;
end

  gNode(:,1) = S(1).data.X;
  gNode(:,2) = S(1).data.Y;