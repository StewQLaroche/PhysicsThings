function [value] = christoffel(l,m,n)
%r coordinate is 1, theta coordinate is 2

value=0;
syms r w;

metric=[1,0;0,r^2*w^2];
inverse=[1,0;0,(1/r^2*w^2)];

d=[0,0;0,2*r*w^2];
d(:,:,2)=[0,0;0,0];

for s=1:2
    value=value+0.5*inverse(l,s)*(d(m,s,n)+d(n,s,m)-d(s,m,n));
end


end

