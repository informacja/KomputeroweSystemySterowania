% przepl.m
function [F1, F2]=przepl(U)
global invh u0 u2 lamb1_wrf lamb2_wrf wBA Sr2Rp Rp P0; 
u1=u0+u2*U;
rf1=(lamb1_wrf)*(wBA/u1+1);
if(invh)
    rf2=(lamb2_wrf)*(wBA/(1-u1)+1);
    c=sqrt(rf1/rf2);
    fR=c*(1-u1)/u1;
    Rf1=Rp*(Sr2Rp*rf1/u1^2+(1+fR)^2);
else
    Rf1=Rp*(Sr2Rp*rf1/u1^2+1);
    fR=0;
end
F1=sqrt(P0/Rf1); F2=F1*fR;
% 
% fR=sqrt(r1/r2)*(1-u1)/u1;
% R(i)=S^2*r1/u1^2+A*(1+fR);
% F1=wF*sqrt(P0/R); F2=F1*fR;

