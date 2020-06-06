% opor.m
function parProc(Unom)
global S r1 r2 A wF wF1 wF2 P0 K1 K2 Kob;
S=0.2; r1=10; r2=100; 
A=1e-3; P0=0.9; 
Unom=1; dU=0.1
u=[Unom-dU/2 Unom+dU/2];
for(i=1:2)
    u1=0.4+0.3*u(i);
    fR=sqrt(r1/r2)*(1-u1)/u1;
    R1=S^2*r1/u1^2+A*(1+fR); 
    f1(i)=sqrt(P0/R1); f2(i)=f1(i)*fR;
    f(i)=f1(i)+f2(i); 
end
ko=(f(2)-f(1))/dU; k1=(f1(2)-f1(1))/dU; k2=(f2(2)-f2(1))/dU; 
wF=Kob/ko; wF1=K1/k1; wF2=-K2/k2;
