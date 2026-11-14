clc
clear

[No_points,phi_points,lam_points,h_points]=textread('Test_Data.txt','%d%f%f%f');

% semi-major axis and semi-minor axis of the reference ellipsoid GRS80
a=6378137;
f=1/298.257222101;
b=a-a*f;

% Geocentric constant GRS80
GM=3986005*10^8;

% Angular velocity of the Earth GRS80
omega=7.292115*10^(-5);

% normalized spherical harmonic coefficient of second order GRS80
C2=-108263*10^(-8)/sqrt(5);

nmax=2190;

for i=1:length(No_points)
    
    i
    
    phi=phi_points(i);
    lam=lam_points(i);
    h=h_points(i);
    
    [xi_con,eta_con]=Conventional_xi_eta(phi,lam,h,a,b,GM,C2,omega,nmax);
    
    [xi_exa,eta_exa]=Exact_xi_eta(phi,lam,h,a,b,nmax);
           
    Conventional_Exact_xi_eta_Test_Data(i,:)=[No_points(i) phi lam h xi_con eta_con xi_exa eta_exa];
    
end 

save Conventional_Exact_xi_eta_Test_Data Conventional_Exact_xi_eta_Test_Data
