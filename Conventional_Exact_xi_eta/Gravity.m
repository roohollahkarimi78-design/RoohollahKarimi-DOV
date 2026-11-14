% This function computes the actual gravity vector

% Inputs
% phi: Geodetic latitude in degrees -90<=phi<=90
% lam: Geodetic longitude in degrees 0<=lam<360
% h: Geodetic height in meters
% a: Semi-major axis of the reference ellipsoid in meters, for GRS80 a=6378137
% b: Semi-minor axis of the reference ellipsoid in meters, for GRS80 b=6356752.3141
% nmax: Maximum degree of expansion

% Outputs
% g_phis: Componet of actual gravity vector along spherical latitude
% g_lams: Componet of actual gravity vector along spherical longitude 
% g_r: Componet of actual gravity vector along spherical radius 
% g: Actual gravity

function [g_phis,g_lams,g_r,g]=Gravity(phi,lam,h,a,b,nmax)

[X,Y,Z]=philamh_XYZ(phi,lam,h,a,b);
[phis,lams,r]=XYZ_phislamsr(X,Y,Z);

% Radius of the reference sphere EIGEN6C4
R=6378136.460;        
 
% Geocentric constant
GM=3986004.415*10^8; 
 
% Angular velocity of the Earth
omega=7.292115*10^(-5);
 
% Loading geopotential coefficients of EIGEN6C4 up to degree and order 2190
load Cnm_EIGEN_6C4_zero_tide Cnm_EIGEN_6C4_zero_tide
load Snm_EIGEN_6C4_zero_tide Snm_EIGEN_6C4_zero_tide
 
% Renaming the variables for the geopotential coefficients
Cnm=Cnm_EIGEN_6C4_zero_tide(1:nmax+1,1:nmax+1);
Snm=Snm_EIGEN_6C4_zero_tide(1:nmax+1,1:nmax+1);
 
% Converting phi and lam from degree to radian
phis=phis*pi/180;
lams=lams*pi/180;
 
% Computing Legendre function of the fist kind and its first derivative
[P,P1]=P_P1_nm(pi/2-phis,nmax);
 
% Computing sin(m_lam) and cos(m_lam)
[co,si]=CosSin(lams,nmax);
 
% Defining the variable n
n=0:nmax;
 
for i=1:length(lams)
    
    m(i,:)=n;
    
end
 
% Derivatives of the gravitational potential
dVg_dphis=-(GM/r)*(R/r).^n*(P1.*Cnm*co'+P1.*Snm*si');
dVg_dlams=(GM/r)*(R/r).^n*(-P.*Cnm*(si.*m)'+P.*Snm*(co.*m)');
dVg_dr=-(GM/r^2)*((n+1).*(R/r).^n)*(P.*Cnm*co'+P.*Snm*si');
 
% Derivatives of the centrifugal potential
dVc_dphis=-omega^2*r^2*cos(phis)*sin(phis);
dVc_dlams=0;
dVc_dr=omega^2*r*cos(phis)^2;
 
% Metric tensor components of the spherical coordinates
gphisphis=r^2;
glamslams=r^2*cos(phis)^2;
grr=1;
 
% Gravity acceleration vector
g_phis=(dVg_dphis+dVc_dphis)/sqrt(gphisphis);
g_lams=(dVg_dlams+dVc_dlams)/sqrt(glamslams);
g_r=(dVg_dr+dVc_dr)/sqrt(grr);
 
% Gravity acceleration
g=sqrt(g_phis.^2+g_lams.^2+g_r.^2);
