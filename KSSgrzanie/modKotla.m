function [y,Twg,Twz,qwy,Fw,Fwy,rog,roz]=modKotla(T,Twg,Twz,qkal)
global wExp ro0 cw Kkal Fvwz;
%CTwz=(qkal/(Fwz*cw)+Twz)/(1+wExp*(Twz-20));
%Twg0=CTwz*(1-20*wExp)/(1-CTwz*wExp);
%rog0=ro0/(1+wExp*(Twg-20));
% (Twg-Twz)=(Twg-T)*(1-exp(-Kkal/(Fw*cw)))]
if(Twg<=T) Twg=T+1.e-6; end
if(Twz>=Twg) Twz=(Twg+T)/2; end
rog=ro0/(1+wExp*(Twg-20));
roz=ro0/(1+wExp*(Twz-20));
% Fw*cw*(Twg-Twz)=Fv*cw*(Twg*rog-Twz*roz);
Fw=Fvwz*(rog*Twg-Twz*roz)/(Twg-Twz);
qK=exp(-Kkal/(Fw*cw));
DTgz=(Twg-T)*(1-qK);
Twz=Twg-DTgz;
roz1=ro0/(1+wExp*(Twz-20));
Fwy=Fvwz*(rog*Twg-Twz*roz)/DTgz;
qK=exp(-Kkal/(Fwy*cw));
qwy=Fwy*cw*(Twg-T)*(1-qK);
y=qkal-qwy;
return;
