function [] = bonding(a)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


a=0.529;
n=1;

for R=0:0.1:5
    
    x=R/a;
    
    X=(1+x)*exp(-x);
    D=(1/x)-((1+(1/x))*exp(-2*x));
    I=exp(-x)*(1+x+(1/3)*x^2);

    bonding_F(n)=1+2*((D+X)/(1+I))-(2/x);
    antibonding_F(n)=1+2*((D-X)/(1-I))-(2./x);

    %bonding(n)=-13.6*(1+bonding_F(n));
    %antibonding(n)=13.6*(1+antibonding_F(n));
    n=n+1;
    
end

hold on
plot(0:0.1:5,(-bonding_F+1)*13.6);
plot(0:0.1:5,(-antibonding_F+1)*13.6);
plot(0:0.1:5,0,'r-');
xlabel('Proton separation in Angstroms');
ylabel('Total energy of system in eV');
axis([0.25 4.5 -2 4]);

end

