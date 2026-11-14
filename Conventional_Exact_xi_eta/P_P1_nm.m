% This function computes the Legendre functions of the first kind and its
% first derivative

% Inputs
% t: Spherical co-latitude in radians, t=pi/2-phis, 0<=t<=pi
% nmax: Maximum degree of expansion

% Outputs
% P: Legendre functions of the first kind
% P1: First derivative of the Legendre functions of the first kind

function [P,P1] = P_P1_nm(t,nmax)

P(nmax+1,nmax+1)=0;
P1(nmax+1,nmax+1)=0;

f11=sqrt(3);
f12=sqrt(5);
f21=sqrt(15);
f22=sqrt(15)/2;

P(1,1)=1;
P(2,1)=f11*cos(t);
P(2,2)=f11*sin(t);
P(3,1)=f12*(1.5*(cos(t))^2-0.5);
P(3,2)=f21*cos(t)*sin(t);
P(3,3)=f22*(sin(t))^2;

for i=3:nmax
    fn1=sqrt(2*i+1)/sqrt(2*i);
    P(i+1,i+1)=fn1*sin(t)*P(i,i);    
    fn2=sqrt(2*i+1)/sqrt(2*(i-1));
    P(i+1,i)=fn2*sin(t)*P(i,i-1);    
    for j=0:i-2	        
        fnm1=sqrt(4*i^2-1)/sqrt(i^2-j^2);
        fnm2=sqrt((2*i+1)*((i-1)^2-j^2))/sqrt((i^2-j^2)*(2*i-3));
        P(i+1,j+1)=fnm1*cos(t)*P(i,j+1)-fnm2*P(i-1,j+1);        
    end    
end

P1(1,1)=0;

for n=1:nmax    
    P1(n+1,1)=-sqrt(n*(n+1)/2)*P(n+1,2);    
end

P1(2,2)=P(2,1);

for n=2:nmax    
    P1(n+1,n+1)=sqrt(n*(n+1)-n*(n-1))/2*P(n+1,n);    
end

for n=2:nmax
    P1(n+1,2)=sqrt(n*(n+1)/2)*P(n+1,1)-sqrt(n*(n+1)-2)/2*P(n+1,3);
    for m=2:n-1
        P1(n+1,m+1)=sqrt(n*(n+1)-m*(m-1))/2*P(n+1,m)-sqrt(n*(n+1)-m*(m+1))/2*P(n+1,m+2);        
    end
end