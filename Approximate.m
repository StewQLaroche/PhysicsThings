function [gravitational_net] = Approximate(in)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
out=in;

gravitational_in=[0];
gravitational_in1=[0];
gravitational_in2=[0];

scale=0;
radial_scale=[0.1:0.1:25.5];
value=0;
G=4.50*10^-39; %in kpc^3/m_solar*s^2
kpc_to_m=3.09*10^19; %converts kiloparsecs to meters
Lnaught=2.25*10^8; %in m_solar/kpc^2
Dnaught=2.25*10^8;
k=0.2; %in 1/kpc
counter=0;
hold=0;

if in==1;
    %computes inner disk with local area cut out
    for rnaught=[0.1:0.1:25.5];
        for theta=0.1:0.1:pi;
            for r=0:0.1:rnaught;
                scale=(rnaught-r*cos(theta));
                scale=scale/(((r^2+rnaught^2-2*r*rnaught*cos(theta)))^0.5);
                hold=(Lnaught*G*scale*r*exp(-k*r));
                hold=hold/(r^2+rnaught^2-2*r*rnaught*cos(theta));
                hold=hold*0.1*0.1;
                value = hold;
            end
        end
        counter=counter+1;
        gravitational_in1(counter)=value;
        value=0;
    end

    counter=0;

    for rnaught=[0.1:0.1:25.5];
        for theta=0:0.1:0.1;
            for r=0:0.1:(rnaught-0.1);
                scale=(rnaught-r*cos(theta));
                scale=scale/(((r^2+rnaught^2-2*r*rnaught*cos(theta)))^0.5);
                hold=(Lnaught*G*scale*r*exp(-k*r));
                hold=hold/(r^2+rnaught^2-2*r*rnaught*cos(theta));
                hold=hold*0.1*0.1;
                value = hold;
            end
        end
        counter=counter+1;
        gravitational_in2(counter)=value;
        value=0;
    end

    gravitational_in=2*(gravitational_in1+gravitational_in2)*kpc_to_m;
    counter=0;

    %computes now with added dark matter
    L=length(gravitational_in);
    gravitational_dmd1=zeros(1,L);
    gravitational_dmd2=zeros(1,L);
    
    %(Dnaught*exp(((r-10)^2)/(-18))*exp(r*0.1))
    %(Dnaught*exp(r*0.07)-Dnaught+Dnaught*exp(((r-10)^2)/(-30)))
 
    for rnaught=[0.1:0.1:25.5];
        for theta=0.1:0.1:pi;
            for r=0:0.1:rnaught;
                scale=(rnaught-r*cos(theta));
                scale=scale/(((r^2+rnaught^2-2*r*rnaught*cos(theta)))^0.5);
                hold=((Lnaught*exp(-k*r)+(Dnaught*exp(r*0.07)-Dnaught+Dnaught*exp(((r-10)^2)/(-30))))*G*scale*r);
                hold=hold/(r^2+rnaught^2-2*r*rnaught*cos(theta));
                hold=hold*0.1*0.1;
                value = hold;
            end
        end
        counter=counter+1;
        gravitational_dmd1(counter)=value;
        value=0;
    end

    counter=0;

    for rnaught=[0.1:0.1:25.5];
        for theta=0:0.1:0.1;
            for r=0:0.1:(rnaught-0.1);
                scale=(rnaught-r*cos(theta));
                scale=scale/(((r^2+rnaught^2-2*r*rnaught*cos(theta)))^0.5);
                hold=((Lnaught*exp(-k*r)+(Dnaught*exp(r*0.07)-Dnaught+Dnaught*exp(((r-10)^2)/(-30))))*G*scale*r);
                hold=hold/(r^2+rnaught^2-2*r*rnaught*cos(theta));
                hold=hold*0.1*0.1;
                value = hold;
            end
        end
        counter=counter+1;
        gravitational_dmd2(counter)=value;
        value=0;
    end

    gravitational_dmd=2*(gravitational_dmd1+gravitational_dmd2)*kpc_to_m;
    counter=0;

    gravitational_net=gravitational_in;
    
    subplot(3,2,2), plot([0.1:0.1:25.5],gravitational_dmd);
    subplot(3,2,2), xlabel('Kiloparsecs');
    subplot(3,2,2), ylabel('Gravitational Acceleration with DMD in m/s^2');
    
    subplot(3,2,3), plot([0.1:0.1:25.5],gravitational_net);
    subplot(3,2,3), xlabel('Kiloparsecs');
    subplot(3,2,3), ylabel('Net Gravitational Acceleration in m/s^2');
    
    velocity=(gravitational_net.*radial_scale.*kpc_to_m).^0.5;
    velocity=velocity/1000; %now in km/s
    subplot(3,2,4), plot([0.1:0.1:25.5],velocity);
    subplot(3,2,4), xlabel('Kiloparsecs');
    subplot(3,2,4), ylabel('Expected Radial Velocity in km/s');
    
    velocity_dmd=(gravitational_dmd.*radial_scale.*kpc_to_m).^0.5;
    velocity_dmd=velocity_dmd/1000; %now in km/s
    subplot(3,2,1), plot([0.1:0.1:25.5],(Dnaught*exp(radial_scale.*0.07)-Dnaught+Dnaught.*exp(((radial_scale-10).^2)/(-30))));
    subplot(3,2,1), xlabel('Kiloparsecs');
    subplot(3,2,1), ylabel('Dark Matter Density in M_s_o_l/kpc^2');
    
end

velocity_observed=velocity;
L=length(velocity_observed);
max_and_location=[0,0];

for i=1:L;
    if max_and_location(2) < velocity_observed(i)
       max_and_location=[i,velocity_observed(i)]; 
    end
end

for i=max_and_location(1):L;
    velocity_observed(i)=(max_and_location(2)+0.02*i);
end

%velocity_observed(max_and_location(1):L)=max_and_location(2);

subplot(3,2,5), plot([0.1:0.1:25.5],velocity_observed);
subplot(3,2,5), xlabel('Kiloparsecs');
subplot(3,2,5), ylabel('Measured Radial Velocity in km/s');

gravitational_needed=((velocity_observed*1000).^2)./(radial_scale*kpc_to_m);
subplot(3,2,6), plot([0.1:0.1:25.5],gravitational_needed);
subplot(3,2,6), xlabel('Kiloparsecs');
subplot(3,2,6), ylabel('Gravitational Acceleration Needed in m/s^2');

%{
dmd=zeros(1,L);
current=0;

for i=1:L;
    while (current < (gravitational_needed(i)/kpc_to_m));
        %true='true'
        
        for theta=0.1:0.1:pi;
            for r=0:0.1:(i*0.1);
                scale=(rnaught-r*cos(theta));
                scale=scale/(((r^2+rnaught^2-2*r*rnaught*cos(theta)))^0.5);
                hold=((Lnaught*exp(-k*r)+dmd(i))*2*G*scale*r);
                hold=hold/(r^2+rnaught^2-2*r*rnaught*cos(theta));
                hold=hold*0.1*0.1;
                current = hold;
            end
        end
        
        for theta=0:0.1:0.1;
            for r=0:0.1:((i-1)*0.1);
                scale=(rnaught-r*cos(theta));
                scale=scale/(((r^2+rnaught^2-2*r*rnaught*cos(theta)))^0.5);
                hold=((Lnaught*exp(-k*r)+dmd(i))*2*G*scale*r);
                hold=hold/(r^2+rnaught^2-2*r*rnaught*cos(theta));
                hold=hold*0.1*0.1;
                current = current + hold;
            end
        end
        
        if (current < (gravitational_needed(i)/kpc_to_m));
            %true2='true2'
            dmd(i) = dmd(i) + 1000000;
        end
        
    end
    %test=dmd(i);
    %test;
end


subplot(3,2,2), plot([0.1:0.1:25.5],dmd);
subplot(3,2,2), xlabel('Kiloparsecs');
subplot(3,2,2), ylabel('Dark Matter Density Distribution in M_sol/kpc^2');

%}






end

