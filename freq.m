function [freq,f_err] = freq(x,y,z,x_err,y_err,z_err,n_x,n_y,n_z,temp,temp_err)
%Given dims. (in cm) and desired mode config, returns needed freq and err.

x=x/100;
y=y/100;
z=z/100;
x_err=x_err/100;
y_err=y_err/100;
z_err=z_err/100;

T=temp;
C=[1.40238744*10^3,5.03836171,5.81172916*10^-2,3.34638117*10^-4,1.48259672*10^-6,3.16585020*10^-9];
v=(C(1))+(C(2)*T)-((C(3))*T^2)+((C(4))*T^3)-((C(5))*T^4)+((C(6))*T^5);
v_err=(C(2)-2*C(3)*T+3*C(4)*T^2-4*C(5)*T^3+5*C(6)*T^4)*temp_err;

A=sqrt((n_x*pi/x)^2+(n_y*pi/y)^2+(n_z*pi/z)^2);

freq=(v/(2*pi))*A;

f_err=0;
holder=(A*v_err/(2*pi))^2;
holder=holder+(n_x*v*x_err/(2*A*x^3))^2;
holder=holder+(n_y*v*y_err/(2*A*y^3))^2;
holder=holder+(n_z*v*z_err/(2*A*z^3))^2;
f_err=sqrt(holder);

end

