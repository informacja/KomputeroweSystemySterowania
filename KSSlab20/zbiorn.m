% przepl.m
%function [F1, F2]=zbiorn(u1)
global WspXn rzad rz1 Tob invh u0 u2 lamb1_wrf lamb2_wrf wBA Sr2Rp Rp P0 Sr A B K1 K2 Kob; 
clear U1 R F1 F2; 
Nielin=1; rzad=5; invh=1; Unom=1; 
Kob=2; ToSum=1; Tob=ToSum/rzad; wKT=2.5; 
wK=wKT*Tob/Tob; %wKT=K1*T2/(K2*T1);
K2=Kob/(wK-1); K1=wK*K2;
[F1nom,F2nom,Fnom,X0]=parProc(Unom,Kob,wKT,ToSum,Tob,Tob,Nielin,invh); 
% ................................
dU=0.01;
U1=[0:dU:Unom 1+dU:dU:1.8]; %U1=[0.4:0.01:0.9]; 
n1=find(U1==Unom); 
lU=length(U1);
for(i=1:lU) 
    [F1(i),F2(i)]=przepl(U1(i)); F(i)=F1(i)+F2(i); 
end
Km=diff(F)/dU; Km1=diff(F1)/dU; Km2=-diff(F2)/dU;
Km=[Km(1) Km]; Km1=[Km1(1) Km1]; Km2=[Km2(1) Km2];
figure(323); 
subplot(2,1,1); 
plot(U1,F,'k',U1(n1),F(n1),'ko');
if(invh)
    hold on;
    plot(U1,F1,'b',U1,F2,'r',U1(n1),F1(n1),'bo',U1(n1),F2(n1),'ro'); hold off;
    txl=sprintf('Przeplywy F(k), F1(b), F2(r); U_{nom}=%.2f F_{nom}=%.2f[m^3/h]=%.3f[l/s] F_{1nom}=%.2f[l/s] F_{2nom}=%.2f[l/s] v=%.1f[cm/s]',Unom, F(n1),F(n1)/3.6,F1(n1)/3.6,F2(n1)/3.6,F(n1)/Sr/36);    
else txl=sprintf('Przeplyw F[m^3/f]; U_{nom}=%.2f F_{nom}=%.2f[m^3/h]=%.3f[l/s] v=%.1f[cm/s]',Unom, F(n1),F(n1)/3.6,F(n1)/Sr/36);    
end
xlabel(txl); axis('tight'); 
if(invh==0) 
    subplot(2,1,2); txl=sprintf('Wzmocnienie K: K_{nom}=%.2f',Km(n1)); 
    plot(U1,Km,'k',U1(n1),Km(n1),'ko'); 
else
    subplot(2,2,3); txl=sprintf('Wzmocnienia K(k), K_1(b) K_2(r): K_{nom}=%.2f K_{1nom}=%.2f K_{2nom}=%.2f',Km(n1),Km1(n1),Km2(n1)); 
    plot(U1,Km1,'b',U1,Km2,'r',U1,Km,'k',...
         U1(n1),Km1(n1),'bo',U1(n1),Km2(n1),'ro',U1(n1),Km(n1),'ko'); axis('tight'); 
end
axis('tight'); xlabel(txl);   
if(invh)
    subplot(2,2,4); 
    plot(U1,Km1./Km2,'k',U1(n1),Km1(n1)/Km2(n1),'ko'); axis('tight'); 
    xlabel(sprintf('Stosunek K_1/K_2; K_{1nom}/K_{2nom}=%.3f',Km1(n1)/Km2(n1)));   
end
