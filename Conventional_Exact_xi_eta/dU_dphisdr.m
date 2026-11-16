% This function computes derivatives of the normal gravity potential with
% respect to spherical latitude and spherical radius

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
% dU_dphis: Derivative of normal gravity potential with respect to spherical latitude
% dU_dr: Derivative of normal gravity potential with respect to
% spherical radius

function [dU_dphis,dU_dr]=dU_dphisdr(phi,lam,h,a,b,GM,C2,omega,nmax)

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

% Derivatives of the normal gravity potential
dU_dphis=dVg_dphis+dVc_dphis;
dU_dr=dVg_dr+dVc_dr;
