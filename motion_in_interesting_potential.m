function [] = motion_in_interesting_potential(a,b,x_0,y_0,dx_0,dy_0)
%tracks point mass in a potential of interest

z=@(x,y) a*(sin((x-pi/2)+y)+sin((x-pi/2)-y))+(x.^2+y.^2)*b;
%[dz_x,dz_y]=[2*(b*x+a*cos(y)*sin(x)),2*(b*y+a*cos(x)*sin(y))];




end

