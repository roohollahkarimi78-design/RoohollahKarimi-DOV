% This function computes the zonal harmonic coefficients of even degree of the normal
% gravity field (see Eq. (2-170) in Hofmann and Mortiz 2005)

% Inputs
% a: Semi-major axis of the reference ellipsoid in meters, for GRS80 a=6378137
% b: Semi-minor axis of the reference ellipsoid in meters, for GRS80 b=6356752.3141
% C2: Normalized zonal harmonic coefficient of the second degree of the
% reference ellipsoid, for GRS80 C2=-108263*10^(-8)/sqrt(5)
% nmax: Maximum degree of expansion

% Outputs
% C: Zonal harmonic coefficients of even degree of the normal gravity field

function C=C2n(a,b,C2,nmax)

E=sqrt(a^2-b^2);

e=E/a;

C2un=C2*sqrt(5);

C(1,1)=1;

C(2,1)=C2un/sqrt(5);

for i=2:nmax/2
    
    Cun=(-1)^i*3*e^(2*i)/((2*i+1)*(2*i+3))*(1-i-5*i*C2un/(e^2));
    
    C(i+1,1)=Cun/sqrt(4*i+1);
    
end
