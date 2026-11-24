% Ths function computes sine and cosine functions

% Inputs
% L: Spherical longitude in radians 0<=L<2*pi
% nmax: Maximum degree of expansion

%Outputs
% co: Cosine function cos(m*L)
% si: Sine function sin(m*L)

function [co,si] = CosSin(L,nmax)

LL=length(L);
 
si(1:LL,1:nmax+1)=0;
co(1:LL,1:nmax+1)=0;
 
si(1:LL,1)=0;
co(1:LL,1)=1;
 
si(1:LL,2)=sin(L);
co(1:LL,2)=cos(L);
 
si1=sin(L);
co1=cos(L);
 
si2=2*cos(L).*sin(L);
co2=2*cos(L).^2-1;
 
si(1:LL,3)=si2;
co(1:LL,3)=co2;
 
for i=3:nmax
    co3=2*cos(L).*co2-co1;
    si3=2*cos(L).*si2-si1;
    si1=si2;
    si2=si3;
    co1=co2;
    co2=co3;    
    si(1:LL,i+1)=si2;
    co(1:LL,i+1)=co2;    
end
