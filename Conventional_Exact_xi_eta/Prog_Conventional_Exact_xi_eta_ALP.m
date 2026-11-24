clc
clear

load ALP ALP

No_points=ALP(:,1);
phi_points=ALP(:,2);
lam_points=ALP(:,3);
h_points=ALP(:,4);

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
    
    No=No_points(i);
    phi=phi_points(i);
    lam=lam_points(i);
    h=h_points(i);
    
    [xi_con,eta_con]=Conventional_xi_eta(phi,lam,h,a,b,GM,C2,omega,nmax);
    
    [xi_exa,eta_exa]=Exact_xi_eta(phi,lam,h,a,b,nmax);
           
    Conventional_Exact_xi_eta_ALP(i,:)=[No phi lam h xi_con eta_con xi_exa eta_exa];
    
end 

min_xi_con=min(Conventional_Exact_xi_eta_ALP(:,5))
max_xi_con=max(Conventional_Exact_xi_eta_ALP(:,5))
mean_xi_con=mean(Conventional_Exact_xi_eta_ALP(:,5))
std_xi_con=std(Conventional_Exact_xi_eta_ALP(:,5))
rms_xi_con=sqrt(mean(Conventional_Exact_xi_eta_ALP(:,5).^2))

min_eta_con=min(Conventional_Exact_xi_eta_ALP(:,6))
max_eta_con=max(Conventional_Exact_xi_eta_ALP(:,6))
mean_eta_con=mean(Conventional_Exact_xi_eta_ALP(:,6))
std_eta_con=std(Conventional_Exact_xi_eta_ALP(:,6))
rms_eta_con=sqrt(mean(Conventional_Exact_xi_eta_ALP(:,6).^2))

min_xi_exa=min(Conventional_Exact_xi_eta_ALP(:,7))
max_xi_exa=max(Conventional_Exact_xi_eta_ALP(:,7))
mean_xi_exa=mean(Conventional_Exact_xi_eta_ALP(:,7))
std_xi_exa=std(Conventional_Exact_xi_eta_ALP(:,7))
rms_xi_exa=sqrt(mean(Conventional_Exact_xi_eta_ALP(:,7).^2))

min_eta_exa=min(Conventional_Exact_xi_eta_ALP(:,8))
max_eta_exa=max(Conventional_Exact_xi_eta_ALP(:,8))
mean_eta_exa=mean(Conventional_Exact_xi_eta_ALP(:,8))
std_eta_exa=std(Conventional_Exact_xi_eta_ALP(:,8))
rms_eta_exa=sqrt(mean(Conventional_Exact_xi_eta_ALP(:,8).^2))

min_dxi=min(Conventional_Exact_xi_eta_ALP(:,7)-Conventional_Exact_xi_eta_ALP(:,5))*1000
max_dxi=max(Conventional_Exact_xi_eta_ALP(:,7)-Conventional_Exact_xi_eta_ALP(:,5))*1000
mean_dxi=mean(Conventional_Exact_xi_eta_ALP(:,7)-Conventional_Exact_xi_eta_ALP(:,5))*1000
std_dxi=std(Conventional_Exact_xi_eta_ALP(:,7)-Conventional_Exact_xi_eta_ALP(:,5))*1000
rms_dxi=sqrt(mean((Conventional_Exact_xi_eta_ALP(:,7)-Conventional_Exact_xi_eta_ALP(:,5)).^2))*1000

min_deta=min(Conventional_Exact_xi_eta_ALP(:,8)-Conventional_Exact_xi_eta_ALP(:,6))*1000
max_deta=max(Conventional_Exact_xi_eta_ALP(:,8)-Conventional_Exact_xi_eta_ALP(:,6))*1000
mean_deta=mean(Conventional_Exact_xi_eta_ALP(:,8)-Conventional_Exact_xi_eta_ALP(:,6))*1000
std_deta=std(Conventional_Exact_xi_eta_ALP(:,8)-Conventional_Exact_xi_eta_ALP(:,6))*1000
rms_deta=sqrt(mean((Conventional_Exact_xi_eta_ALP(:,8)-Conventional_Exact_xi_eta_ALP(:,6)).^2))*1000

save Conventional_Exact_xi_eta_ALP Conventional_Exact_xi_eta_ALP
