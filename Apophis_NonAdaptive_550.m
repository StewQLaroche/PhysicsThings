function Apophis_NonAdaptive_550

%Simplified and inaccurate code to calculate the orbit of the astroid
%aphosis using Euler's method with nonadpative time stepping.

rr=zeros(400000,1);
plot_every = 60*60*24*8;
time_step=3600;
time_step_days=time_step/(24*60*60);
time_step_years=time_step/(24*60*60*365);

startdate=734068;
gravitational_constant=6.672e-11;
mass_earth=5.98e24;
mass_sun=1.991e30;
mass_apophis=2.7e10;
mass_venus=4.88e24;

radius_e=1.496e11;

y=zeros(12,1);
y(1)=1.4960e11;       %x position of the Earth at the start_date given in m
y(2)=0.0;             %y position of the Earth at the start_date given in m
y(3)=0.0;             %x velocity of the Earth at the start_date given in m/s
y(4)=2.9783e4;        %y velocity of the Earth at the start_date given in m/s
y(5)=-4.128199e10;    %x position of the Venus at the start_date given in m
y(6)=9.980061e10;     %y position of the Venus at the start_date given in m
y(7)=-3.232230e4;     %x velocity of the Venus at the start_date given in m/s
y(8)=-1.336962e4;     %y velocity of the Venus at the start_date given in m/s
y(9)=-1.468739e11;    %x position of the Aphosis at the start_date given in m
y(10)=-2.230613e10;   %y position of the Aphosis at the start_date given in m
y(11)=9.127491e3;     %x velocity of the Aphosis at the start_date given in m/s
y(12)=-2.724475e4;    %y velocity of the Aphosis at the start_date given in m/s

gamma=y;
GM=gravitational_constant*mass_sun;
GE=gravitational_constant*mass_earth;
GV=gravitational_constant*mass_venus;
GA=gravitational_constant*mass_apophis;


%Explicit Euler
figure(1)
k=1;
j=0;
t=0;
pp=0;
rr(1)=1e12;
c_date=startdate;
count=0;

while (t<35)
j=j+time_step;
t=t+time_step_years;
c_date=c_date+time_step_days;
count=count+1;



% The only section of the code you should change is below here
%d=y+time_step*ff(y);
%y=y+(time_step*ff(y))+(ff(d)/2);

k1=ff(y);
k2=ff(y+(time_step/2)*k1);
k3=ff(y+(time_step/2)*k2);
k4=ff(y+time_step*k3);
y=y+((time_step/6)*(k1+2*k2+2*k3+k4));

% The only section of the code you should change is above here

rr(count)=norm([y(1)-y(9),y(2)-y(10)]);

% Everything inside the following if statement is for graphics only
    if ~mod(j,(plot_every))
        plot(y(1)/radius_e,y(2)/radius_e,'k*','linewidth',3);
        hold on;
        plot(y(9)/radius_e,y(10)/radius_e,'k*','linewidth',3);
        plot(y(5)/radius_e,y(6)/radius_e,'k*','linewidth',3);
        hold off;
        axis equal;
        title(datestr(c_date,1));
        %c_date 
        %datestr(c_date,1)
        if j<365*24*60*60+2*plot_every;
            g(1,k)=y(1);
            g(2,k)=y(2);
            g(3,k)=y(9);
            g(4,k)=y(10);
            g(5,k)=y(5);
            g(6,k)=y(6);
            k=k+1;
            hold on;
            plot(0,0,'yo','linewidth',5);
            plot(y(1)/radius_e,y(2)/radius_e,'bo','linewidth',3);
            plot(y(9)/radius_e,y(10)/radius_e,'go','linewidth',3);
            plot(y(5)/radius_e,y(6)/radius_e,'mo','linewidth',3);
            hold off;
        else
            hold on;
            plot(0,0,'yo','linewidth',5);
            plot(g(1,:)/radius_e,g(2,:)/radius_e,'b','linewidth',2);
            plot(g(3,:)/radius_e,g(4,:)/radius_e,'g','linewidth',2);
            plot(g(5,:)/radius_e,g(6,:)/radius_e,'m','linewidth',2);
            plot(y(1)/radius_e,y(2)/radius_e,'k*','linewidth',3);
            plot(y(1)/radius_e,y(2)/radius_e,'bo','linewidth',3);
            plot(y(9)/radius_e,y(10)/radius_e,'go','linewidth',3);
            plot(y(5)/radius_e,y(6)/radius_e,'mo','linewidth',3);
            hold off;
        end    
        axis([-1.5,1.5,-1.5,1.5]);
        drawnow;
        if (round(c_date/100)==7353) pause(0.05); end
        if (round(c_date/200)==round(7382/2)) pause(0.05); end
        if (round(c_date/1000)==round(741)) pause(0.05); end
    end
% Everything inside the above if statement is for graphics only


end


figure(3)
semilogy(rr)

[a,b]=min(rr(1:count));
datestr(b*time_step_days+startdate,1)
a


    function zf=ff(yf)
        zf=0*yf;
        rfes=(yf(1)^2+yf(2)^2)^1.5;
        rfea=(((yf(1)-yf(9))^2)+((yf(2)-yf(10))^2))^1.5;
        rfvs=(yf(5)^2+yf(6)^2)^1.5;
        rfva=(((yf(5)-yf(9))^2)+((yf(6)-yf(10))^2))^1.5;
        rfas=(yf(9)^2+yf(10)^2)^1.5;
        zf(1)=yf(3);
        zf(2)=yf(4);
        zf(5)=yf(7);
        zf(6)=yf(8);
        zf(9)=yf(11);
        zf(10)=yf(12);
        zf(3)=-yf(1)*GM/(rfes)-(yf(1)-yf(9))*GA/(rfea);
        zf(4)=-yf(2)*GM/(rfes)-(yf(2)-yf(10))*GA/(rfea);
        zf(7)=-yf(5)*GM/(rfvs)-(yf(5)-yf(9))*GA/(rfva);
        zf(8)=-yf(6)*GM/(rfvs)-(yf(6)-yf(10))*GA/(rfva);
        zf(11)=-yf(9)*GM/(rfas)+(yf(1)-yf(9))*GE/(rfea)+(yf(5)-yf(9))*GV/(rfva);
        zf(12)=-yf(10)*GM/(rfas)+(yf(2)-yf(10))*GE/(rfea)+(yf(6)-yf(10))*GV/(rfva);
    end
end