%function [y T Ts Tok Twz Twg Fw qKal rog roz]=Ogrzew(T,Ts,Twg,Twz,Tzewn,qK)
function [y T Ts Twz Twg Fw qKal rog roz]=Ogrzew(T,Ts,Twz,Twg,Twgt,Tzewn,qK,Fp)
% Twg=X(4:nTwgl); U=Twgt;
global dt Tob KTTs KTsT KTsTzw Kob Taus KTq0 wAqTqT Aq0 AqTg AqT Mpowcw Mpow KqScian KqSciIz MSccw; 
global nTwg nTwg1 Fvwz Fwconst lX; 
Twg(1:nTwg)=Twg(2:nTwg1); Twg(nTwg1)=Twgt; 
if(Fwconst)
    [Twgwy,Twz,qKal,Fw,Fwyx,brak,rog,roz]=KociolFwz(T,Twg(1),Twz,qK);
else
    [Fw,Twz,qKal,Twgwy, brak,rog,roz]=Kociol(T,Twg(1),Twz,qK); % Opis w WymiennikiTeoria.m
end
% ........ Dynamika grzania pomieszcenia T(t):  % Opis w PrzenikanieTeoria.m
% d(Vp*rop*cwp*T)/dt=Fp*cwp*rop*Tzew-Fp*cwp*rop*T+qKal-(T-Ts)*KqScian
% dT/dt=Fp/Vp*Tzew-Fp/Vp*T+qKal/(cwp*rop*Vp)-(T-Ts)*S*lamb/(cwp*rop)/(Vp*dl)
% dT/dt=-Fp/Vp*T+Fp/Vp*Tzew+qKal/Mpowcw-(T-Ts)KqScian/Mpowcw

%dT_dt=-T/Tauw+Ts/Tauw+qKal/Mpowcw;  KTq=1./KqScian; 
%dT_dt=-T/Tauw+Ts/Tauw+qKal*(KqScian/Mpowcw)/KqScian; 

%dT_dt=(-T+Ts+qKal*KTq)/Tauw; 
%qKal=(Aq0+AqTg*Twg+AqT*T)/(1+wAqTqT*T); 
%qKal=+AqT*T/(1+wAqTqT*T)+(Aq0+AqTg*Twg)/(1+wAqTqT*T); 
% dT/dt=-(Fp/Vp+KqScian-AqT/(1+wAqTqT*T))*T+Ts*KqScian+Fp/Vp*Tzew+qKal/Mpowcw

%dT_dt=(-T+Ts+KTq*(Aq0+AqTg*Twg+AqT*T)/(1+wAqTqT*T))/Tauw; 
%dT_dt=(-T+Ts+(KTq*AqTg)*Twg+(KTq*AqT)*T+KTq*Aq0)/(1+wAqTqT*T))/Tauw; 

%KTTgp=(KTq*AqTg)*3600; KTT=(KTq*AqT)*3600; KT0=KTq*Aq0*3600;
%dT_dt=(-T*(1-KTT)+Ts*(1-KTT)/(1-KTT)+KTTgp*(1-KTT)/(1-KTT)*Twg+KT0*(1-KTT)/(1-KTT))/(1+wAqTqT*T))/Tauw; 
%Tob=Tauw/(1-KTT); KTTs=1/(1-KTT); KTTg=KTTgp/(1-KTT); KTq0=KT0/(1-KTT); 
% Dynamika pomieszczenia:
%qKal=(Aq0+AqTg*Twg+AqT*T)*3600/(1+wAqTqT*T);
dT_dt=-(Fp/Mpow)*(T-Tzewn)-(T-Ts)*KqScian/Mpowcw+qKal/Mpowcw;
%Dynamika grzania �cian
%dTs_dt=((T-Ts)*KqScian-(Ts-Tzewn)*KqSciIz)/MSccw;
dTs_dt=-Ts/Taus+(T*KqScian+Tzewn*KqSciIz)/MSccw;
%dT_dt=(-T+Ts*KTTs+Kob*Twg+KTq0)/Tob/(1+wAqTqT*T); 

% Dynamika grzania �cian
%dTs_dt=-Ts/Taus+T*KqScian/MSccw+Tzewn*KqSciIz/MSccw; % Ks/MSccw=1/Taus; KqSciIz/MSccw=(Ks/MSccw)*KqSciIz/Ks=KTsTzw/Taus 
%dTs_dt=-Ts/Taus+T/Taus*KTsT+Tzewn/Taus*KTsTzw; % KTsT=KqScian/Ks; KTsTzw=KqSciIz/Ks 

%[dTs_dt,Ts]=Sciany(Ts,T,Tzewn); %dTs_dt=-Ts+T*KTsT+Tzewn*KTsTzw)/Taus;
%dTs_dt=Sciany(Ts,T,Tzewn); %dTs_dt=-Ts+T*KTsT+Tzewn*KTsTzw)/Taus;
T=T+dT_dt*dt; 
Ts=Ts+dTs_dt*dt; 
y=T; 
