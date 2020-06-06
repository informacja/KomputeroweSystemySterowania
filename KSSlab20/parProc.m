% opor.m
function [F1nom,F2nom,Fnom,X0]=parProc(Unom,Kob,wKT,Tobsum,Tob1,Tob2,Nielin,invh)
%clear all;
global FU0 WspXn rzad rz1 u0 umax u2 lamb1_wrf lamb2_wrf wBA Sr2Rp Rp P0 Sr A B Swy;
global S1 Sf C1 Cf Tob T2; % K1 K2 Kob; 
%Nielin=1; rzad=5; 
%ToSum=1; Tob=ToSum/rzad;
wK=wKT*Tob1/Tob2; Tob=Tob1; T2=Tob; %wKT=K1*T2/(K2*T1);
K2=Kob/(wK-1); K1=wK*K2;
WspXn=Nielin-1; 
if(invh==0) K1=Kob; K2=0; end
X0=ones(rzad,1); 
rz1=rzad-invh; 
if(Nielin==0) 
    F1nom=K1*Unom; F2nom=K2*Unom; Fnom=F1nom-F2nom;
    X0=X0*Fnom;  S1=10; Sf=S1; C1=2; Cf=C1;
    return; 
end
% ................ Prawo Darcy'ego: ..................... 
% dP/dl=r/D*v^2; dla turb. r=const
% D=4*S/O; S- pole; O obwód obmywany
% Tu mamy kana³ prostok¹tny o wysok.B i szerok. A podziel. na a1, a2, gdzie a2=A-a1; 
% Opor hydrauliczny takich kanalow wyra¿aj¹ wzory: 
%     rf1=(lamb1*((B/A)/u1+1);
%     rf2=(lamb2*((B/A)/(1-u1)+1);
% Mamy równania:
%  DP=(v1^2)*rf1=(F1^2)*(rf1/(s1)^2)
%  DP=(v2^2)*rf2=(F2^2)*(rf2/(s2)^2); s2=S-s1; s2=S(1-u1); s1=S*u1
%  
%  F2=F1*sqrt(rf1/rf2)*(s2/s1)=F1*sqrt(rf1/rf2)*((1-u1)/u1) 
% Równanie pe³nej hydrauliki ma postaæ:
%  DP=P0-Rp*(F1+F2)^2, st¹d
% 
%  P0-Rp*(F1+F2)^2=F1^2*((rf1/((S^2)*u1^2),
%  P0-(F1^2)*Rp*(1+F2/F1)^2=F1^2*((rf1/(S^2))/u1^2),
%  F2/F1=fr=sqrt(rf1/rf2)*((1-u1)/u1); 
%  P0=(F1^2)*(Rp*((1+fr)^2)+((rf1/(S^2))/u1^2)
%  P0=(F1^2)*Rf; 
%  Rf=Rp*((1+fr)^2+((rf1/(Rp*S^2))/u1^2)
%  F1=sqrt(P0/Rf)
umax=0.92; u0=0.2; u2=(umax-u0)/2;  
% ......................
Sr=0.20e-2; r1=10; r2=100; 
r1=15.9; % to daje Km1/Km2=K1/K2
Rp=1.e-3; P0=0.9*3600^2; 
Unom=1;
% ..... Param. dodatk. .....
%clear A B wrf;
A=0.08; %A=0.09; %A=1.2;
if(invh==0) Sr=Sr/2; A=A/2; end
B=Sr/A; wBA=B/A; 
Sr2Rp=1/((Sr^2)*Rp); 
%wrf=2*wBA^3*A^5; %wrf=(2*B); 
% obliczenie wsp.rf1 i rf2: dla Unom rf1=r1 i rf2=r2
u1=u0+u2*Unom;
lamb1_wrf=r1/(wBA/u1+1); %lamb1=lamb1_wrf*wrf;
lamb2_wrf=r2/(wBA/(1-u1)+1); %lamb2=lamb2_wrf*wrf; 
% ................................
T1=Tob1; T2=Tob1; Tob=T1; 
%WKo=wKT*T1/T2; % to ma byc sta³e
% ...... obliczamy r1 tak, aby wKmTm=KwKT ......
clear U u f1 f2 rr1 K1K2 DwK;
if(invh)
    rr1=[10:0.001:30]; lR=length(rr1); %lR=0;
    dU=0.01;
    U=[Unom-dU Unom];
    u=u0+u2*U; fr=(1-u)./u;
    for(j=1:lR)
        rx=rr1(j);  u1=u(2);
        lamb1_wrf=rx/(wBA/u1+1); %lamb1=lamb1_wrf*wrf;
        lamb2_wrf=r2/(wBA/(1-u1)+1); %lamb2=lamb2_wrf*wrf;
        for(i=1:2)
            u1=u(i);
            rf1=(lamb1_wrf)*(wBA/u1+1);
            rf2=(lamb2_wrf)*(wBA/(1-u1)+1);
            c=sqrt(rf1/rf2);
            fr=c*(1-u1)/u1;
            R1=Rp*(Sr2Rp*rf1/u1^2+(1+fr)^2);
            f1(i)=sqrt(P0/R1); f2(i)=f1(i)*fr;
        end
        Km1=(f1(2)-f1(1))/dU; Km2=-(f2(2)-f2(1))/dU;
        K1K2(j)=Km1/Km2; xr1(j)=rx; yrf(j)=rf1/rf2;
        DwK(j)=(Km1*T2/(Km2*T1))/wKT-1;
        if(DwK(j)>0) r1=rx; co=c; end; %break; end
    end
    u1=u0+u2*Unom;
    lamb1_wrf=r1/(wBA/u1+1); %lamb1=lamb1_wrf*wrf;
    lamb2_wrf=r2/(wBA/(1-u1)+1); %lamb2=lamb2_wrf*wrf;
end
modAB=0;
% .....................................
clear U f1 f2 f F1 F2 F; 
% test pe³nej hydrauliki
du=0.01; 
%U1=[0:du:Unom Unom+du:du:1.8]; %U1=[0.4:0.01:0.9]; 
U1=[Unom Unom+du]; %U1=[0.4:0.01:0.9]; 
n1=find(U1==Unom); 
lU=length(U1);
%cr=sqrt(r1/r2);
for(k=1:2)
    for(i=1:lU) 
        [f1(i),f2(i)]=przepl(U1(i)); f(i)=f1(i)+f2(i); 
    end
    km=diff(f)/du; 
    ko=km(n1);
    wF=Kob/ko; 
    if(k==1)
        Rp0=Rp; Sr0=Sr; A0=A; B0=B; r10=r1; r20=r2;
        Rp=Rp/wF^2; 
        wFx=(wF^(1./4)); 
        if(modAB)
            A=A0*wFx; B=B0*wFx; r1=r10/wF; r2=r20/wF; 
        else wF2=wF^2; r1=r10/wF2; r2=r20/wF2;  
        end
        Sr=A*B; Sr2Rp=1/((Sr^2)*Rp);
        wBA=B/A; 
        %lam1=2*B0*r1; lam2=2*B0*r2; 
        %r1=lam1/(2*B); r2=lam2/(2*B); 
        %wrf=2*wBA^3*A^5;
        % obliczenie wsp.rf1 i rf2: dla Unom rf1=r1 i rf2=r2
        u1=u0+u2*Unom;
        lamb1_wrf=r1/(wBA/u1+1); %lamb1=lamb1_wrf*wrf;
        lamb2_wrf=r2/(wBA/(1-u1)+1); %lamb2=lamb2_wrf*wrf;
    end
end
Fnom=f(n1); F1nom=f1(n1); F2nom=f2(n1); 
% ................................
SzbSum=10*rzad; 
Szb1=SzbSum/rzad; S1=SzbSum/rz1; Sf=S1;
X0(rz1)=1; %(Tsum/rzad)*F1nom/(S1*2);  
for(i=rz1-1:-1:1) X0(i)=WspXn*X0(i+1)+X0(rz1); end
X0(rzad)=X0(rz1); 
C1=2*sqrt(X0(rz1))/T1; %F1nom/S1/sqrt(X0(rz1)); 
S1=F1nom/(C1*sqrt(X0(rz1)));
%C1=c1/S1*sqrt(2*9.81*3600^2);
Swy=C1*S1/sqrt(2*9.81*3600^2); 
Cf=C1; if(invh) Sf=F2nom/(C1*sqrt(X0(rzad))); else Sf=S1; end
a=sqrt(Swy);
[f1(1),f2(1)]=przepl(0); FU0=f1(1)+f2(1); 
return;
% ..................
