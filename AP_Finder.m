function [AP] = AP_Finder(max_level)
%Finds all accidental degenerary pairs
%   Using E=(n_x^2+n_y^2), this program searched for n_x,n_y pairs that
%   have the same energy, while not belonging to the standard 2 two-fold
%   degeneracy conditions. Thus, degenerate states of the nature x,y & y,x
%   are discarded. The program arranges pairs in matrix form, with the
%   first column being energy, and then having energy pairs in other
%   columns with pairs separated with a zero.

AP=[0,0;0,0];
column=1;
for i=1:max_level;          %First, create allowable energies
    for j=i:max_level;
        E=i^2 + j^2;
        AP(E,column)=E;
        for x=1:max_level;
            for y=x:max_level;
                if(x^2 + y^2) == E;     %next, check for pairs that 
                    column=column+1;    %correspond to each energy level
                    AP(E,column)= x;    %the method of indexing prevents
                    column=column+1;    %the normal x,y / y,x degeneracy
                    AP(E,column)= y;    %from appearing
                    column=column+1;
                    AP(E,column)= 0;
                end
            end
        end
        column=1;
    end
end
[m,n]=size(AP);
nothing=zeros(1,n);
for z=m:-1:1;               %this loop gets rids of rows that contain no
    if AP(z,:) == nothing;  %energy level, thus cleaning up the matrix
        AP(z,:)=[];
    end
end
[m,n]=size(AP);
no_pairs=zeros(1,n-3);
for z=m:-1:1;                   %this loop gets rid of energy levels that
    if AP(z,4:n) == no_pairs;   %have no accidental degeneracy, again
        AP(z,:)=[];             %cleaning up the matrix
    end
end
end