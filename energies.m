function [E_z] = energies(n,l,s,Z,m_j)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

j=l+s;

E_0 = (-13.6*Z^2)/n^2;

E_r = (-(E_0^2)/(1.022*10^6))*((4*n/(l+0.5))-3);

E_so = ((E_0^2)/(511*10^3))*(n*(j*(j+1)-l*(l+1)-0.75)/(l*(l+0.5)*(l+1)));

E_fs = ((E_0^2)/(1.022*10^6))*((-4*n/(j+0.5))+3);

E_z = (1+(1/(2*j*(j+1)))*(j*(j+1)-l*(l+1)+0.75))*(5.79*10^-5)*m_j;

E = [E_0, E_r, E_so, E_fs];



end