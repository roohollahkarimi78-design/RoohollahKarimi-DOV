% This function computes the components of the Helmert deflection of the
% vertical using the Conventional method (Jekeli 1999)

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
% xi: north-south component of the Helmert deflection of the vertical in
% arc-seconds
% eta: east-west component of the Helmert deflection of the vertical in
% arc-seconds

function [xi,eta]=Conventional_xi_eta(phi,lam,h,a,b,GM,C2,omega,nmax)

e2=(a^2-b^2)/a^2;

[X,Y,Z]=philamh_XYZ(phi,lam,h,a,b);
[phis,lams,r]=XYZ_phislamsr(X,Y,Z);

[dW_dphis,dW_dlams,dW_dr]=dW_dphisdlamsdr(phi,lam,h,a,b,nmax);
[dU_dphis,dU_dr]=dU_dphisdr(phi,lam,h,a,b,GM,C2,omega,nmax);

[gamma_phis,gamma_r,gamma]=Gravity_Normal(phi,lam,h,a,b,GM,C2,omega,nmax);

xi=-1/(gamma*r)*(dW_dphis-dU_dphis)*180/pi*3600+(1/gamma)*(phi-phis)*(dW_dr-dU_dr)*3600+0.17*(h/1000)*sin(2*phi*pi/180);
eta=-1/(gamma*r*cos(phis*pi/180))*dW_dlams*180/pi*3600;
