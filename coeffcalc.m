function [ y ] = coeffcalc(a )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
a
t=[-1:0.01:1];
    
for i=1:50
    coeffs(i)=(-(-1)^i)/(2*pi*i);
end

for f=1:201
    y(f)=0;
end

for z=1:201
    for j=1:50
    y(z) = y(z)+ coeffs(j)*sin(j*z);
    end
end

%for d=1:2000
 %   y(d)
%end

plot(t,y)
    return

end

