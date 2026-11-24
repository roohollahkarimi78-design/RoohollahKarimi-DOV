% This function computes the components of the Helmert deflection of the
% vertical using the exact method (proposed method)

% Inputs
% phi: Geodetic latitude in degrees -90<=phi<=90
% lam: Geodetic longitude in degrees 0<=lam<360
% h: Geodetic height in meters
% a: Semi-major axis of the reference ellipsoid in meters, for GRS80 a=6378137
% b: Semi-minor axis of the reference ellipsoid in meters, for GRS80 b=6356752.3141
% nmax: Maximum degree of expansion

% Outputs
% xi: north-south component of the Helmert deflection of the vertical in
% arc-seconds
% eta: east-west component of the Helmert deflection of the vertical in
% arc-seconds

function [xi,eta]=Exact_xi_eta(phi,lam,h,a,b,nmax)

[g_phis,g_lams,g_r,g]=Gravity(phi,lam,h,a,b,nmax);

n_Sph=-[g_phis,g_lams,g_r]./g;

[X,Y,Z]=philamh_XYZ(phi,lam,h,a,b);
[phis,lams,r]=XYZ_phislamsr(X,Y,Z);

phi=phi*pi/180;
lam=lam*pi/180;

phis=phis*pi/180;
lams=lams*pi/180;

G=[-sin(phis)*cos(lams) -sin(phis)*sin(lams) cos(phis);
    -sin(lams) cos(lams) 0;
    cos(phis)*cos(lams) cos(phis)*sin(lams) sin(phis)];

n_Car=G\n_Sph';

phiast=asin(n_Car(3));
lamast=2*atan(n_Car(2)/(n_Car(1)+sqrt(n_Car(1)^2+n_Car(2)^2)));

if lamast<0
    lamast=lamast+2*pi;
end

eta=(lamast-lam)*cos(phi);
xi=phiast-phi+(1/2)*(eta^2)*tan(phi);

eta=eta*180/pi*3600;
xi=xi*180/pi*3600;
