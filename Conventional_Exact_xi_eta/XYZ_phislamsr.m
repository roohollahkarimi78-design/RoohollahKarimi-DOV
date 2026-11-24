% This function converts Cartesian coordinates to spherical coordinates

% Inputs
% X: X-component of Cartesian coordinates in meters
% Y: Y-component of Cartesian coordinates in meters
% Z: Z-component of Cartesian coordinates in meters

% Outputs
% phis: Spherical latitude in degrees
% lams: Spherical longitude in degrees
% r: Spherical radius in meters

function [phis,lams,r]=XYZ_phislamsr(X,Y,Z)

% Computing lam_spherical from Cartesian coordinates
lams=atan(abs(Y/X));

if (X>0 & Y>=0)  lams=lams;        end
if (X<0 & Y>0)   lams=pi-lams;     end
if (X<0 & Y<0)   lams=pi+lams;     end
if (X>0 & Y<0)   lams=2*pi-lams;   end
if (X==0 & Y>0)  lams=pi/2;     end
if (X==0 & Y<0)  lams=3*pi/2;   end

% Computing r from Cartesian coordinates
r=sqrt(X^2+Y^2+Z^2);

% Computing phi_spherical from Cartesian coordinates
phis=atan(sqrt(Z^2/(X^2+Y^2)));

if Z<0
    phis=-phis;
end

phis=phis*180/pi;
lams=lams*180/pi;

