% function [dTs_dt Ts]=Sciany(Ts,T,Tzewn)
% Dynamika grzania œcian
%dTs_dt=-Ts/Taus+T*Kqw/Mcs+Tzewn*Kqs/Mcs; % Ks/Mcs=1/Taus; Kqs/Mcs=(Ks/Mcs)*Kqs/Ks=KTsTzw/Taus 
%dTs_dt=-Ts/Taus+T/Taus*KTsT+Tzewn/Taus*KTsTzw; % KTsT=Kqw/Ks; KTsTzw=Kqs/Ks 
function dTs_dt=Sciany(Ts,T,Tzewn) %dTs_dt=-Ts+T*KTsT+Tzewn*KTsTzw)/Taus;
global rzads rzad rps KTsT KTsTzw KTokT KTokTzw KTsTiz KTsTs Taus;
Tp=T; KTsTp=KTsT; 
if(rps==2) % sa okna; rzads - rzad scian z izolacj¹
    % Okna: 1.szy rzad 
    dTs(1)=(-Ts(1)+T*KTokT+Tzewn*KTokTzw)/Taus(1);
end
if(rzads>1)
    Tnast=Ts(rps+1); if(rzads>2) KTsTnast=KTsTs; else KTsTnast=KTsTzw; end
else Tnast=Tzewn; KTsTnast=KTsTzw;
end
dTs_dt(rps:rzad)=zeros(1,rzads); 
%for(r=rps:rzad)
r=rps; 
while(1) %r<=rzad)    
    dTs_dt(r)=(-Ts(r)+Tp*KTsTp+Tnast*KTsTnast)/Taus(r);    
    Tp=Ts(r); r=r+1; if(r>rzad) break; end, 
    if(r==rzad) Tnast=Tzewn; KTsTnast=KTsTzw; KTsTp=KTsTiz; 
    else Tnast=Ts(r+1); if(r==rzad-1) KTsTnast=KTsTiz; else KTsTnast=KTsT; end 
    end
% Izolacja zewnêtrzna: 1.szy rzad; Tiz=Ts(rzad)
% dTiz_dt=(-Tiz+Ts(rzad)*KTsTiz+Tnast*KTizTnast)/Taus(r); 
end 

