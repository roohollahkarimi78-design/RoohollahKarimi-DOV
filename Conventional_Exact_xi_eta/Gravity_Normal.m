% This function computes the normal gravity vector

% Inputs
% phi: Geodetic latitude in degrees -90<=phi<=90
% lam: Geodetic longitude in degrees 0<=lam<360
% h: Geodetic height in meters
% a: Semi-major axis of the reference ellipsoid in meters, for GRS80 a=6378137
% b: Semi-minor axis of the reference ellipsoid in meters, for GRS80 b=6356752.3141
% GM: Geocentric gravitational constant of the reference ellipsoid in m^3s^-2, for GRS80 GM=3986005*10^8
% C2: Normalized zonal harmonic coefficient of the second degree of the
% reference ellipsoid, for GRS80 C2=-108263*10^(-8)/sqrt(5)
% omega: Angular velocity of the Earth in rads^-1, omega=7.292115*10^(-5)
% nmax: Maximum degree of expansion

% Outputs
% gamma_phis: Componet of normal gravity vector along spherical latitude 
% gamma_r: Componet of normal gravity vector along spherical radius 
% gamma: Normal gravity

function [gamma_phis,gamma_r,gamma]=Gravity_Normal(phi,lam,h,a,b,GM,C2,omega,nmax)
 
[X,Y,Z]=philamh_XYZ(phi,lam,h,a,b);
[phis,lams,r]=XYZ_phislamsr(X,Y,Z);

% Normal coefficients 
C=C2n(a,b,C2,nmax);

% Converting phi and lam from degree to radian
phis=phis*pi/180;
 
% Computing Legendre function of the fist kind and its first derivative
[P,P1]=P_P1_nm(pi/2-phis,nmax);
 
% Defining the variable n
n=0:2:nmax;
 
% Derivatives of the normal gravitational potential
dVg_dphis=-(GM/r)*(a/r).^n*(P1(1:2:end,1).*C);
dVg_dr=-(GM/r^2)*((n+1).*(a/r).^n)*(P(1:2:end,1).*C);
 
% Derivatives of the normal centrifugal potential
dVc_dphis=-omega^2*r^2*cos(phis)*sin(phis);
dVc_dr=omega^2*r*cos(phis)^2;
 
% Metric tensor components of the spherical coordinates
gphisphis=r^2;
grr=1;
 
% Normal gravity acceleration vector
gamma_phis=(dVg_dphis+dVc_dphis)/sqrt(gphisphis);
gamma_r=(dVg_dr+dVc_dr)/sqrt(grr);
 
% Normal gravity acceleration
gamma=sqrt(gamma_phis.^2+gamma_r.^2);
