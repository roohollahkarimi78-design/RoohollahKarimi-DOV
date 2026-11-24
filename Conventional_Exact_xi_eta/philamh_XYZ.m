% This function converts geodetic coordinates to Cartesian coordinates

% Inputs
% phi: Geodetic latitude in degrees -90<=phi<=90
% lam: Geodetic longitude in degrees 0<=lam<360
% h: Geodetic height in meters
% a: Semi-major axis of the reference ellipsoid in meters, for GRS80 a=6378137
% b: Semi-minor axis of the reference ellipsoid in meters, for GRS80 b=6356752.3141

% Outputs
% X: X-component of Cartesian coordinates in meters
% Y: Y-component of Cartesian coordinates in meters
% Z: Z-component of Cartesian coordinates in meters

function [X,Y,Z]=philamh_XYZ(phi,lam,h,a,b)

% Second eccentricity
e2=(a^2-b^2)/a^2;

% Geodetic latitude 
phi=phi*pi/180;

% Geodetic longitude
lam=lam*pi/180;

% Normal radius of the curvature
N=a/sqrt(1-e2*(sin(phi))^2);   

% Conversion from geodetic to Cartesian coordinates
X=(N+h)*cos(phi)*cos(lam);       
Y=(N+h)*cos(phi)*sin(lam);       
Z=(N*(1-e2)+h)*sin(phi);  
