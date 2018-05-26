function [A] = zeeman(n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

A=[0 0 0 0 0];
index=1;
for l=0:(n-1)
    for s=-0.5:0.5
        for ml=-l:l
            if(l==0)
                outer=((0.75/3)-1);
            else
                outer=((0.75/3)-((l*(l+1)-ml*s)/(l*(l+0.5)*(l+1))));
            end
            zeeman=(10000*20*5.79*10^-5)*(ml+2*s);
            energy=(13.9/17)*10000*(1/137)^2*outer;
            A(index,:)=[l,ml,s,energy,zeeman];
            index=index+1;
        end
    end
end
A=sortrows(A,5);

end

