function [hPom,Tos,Kos,dels,ds,Tor,Kor,delR,dr,Kob,DtR,nDt,nKpocz]=IdentOb(Ldanych,t,dt,U,Y,Ysk0,tSkok,nSkok,Unom,Upocz,Ypocz,Yt,stdV,DUskok,DYshift,ntr,rzad,Kob,Ts,Tokna)
% Identyfikacja obiektu G(s)=K*exp(-d)/(1+Ts)
% Zadajemy hmin i hmax; 
% =============================================================
global przyr wagi nfig1 alfa delr txPar; % Yr Ys; 
przyr=1; przyr=0; wagi=1; nfig1=1;
figure(nfig1); 
kol='k'; %kol='k.-'; 
subplot(2,1,1); hold off;
%Ysk0=Fnom-Unom*Kob; %Y0lin-FU0;

plot(t,(U-Unom)*Kob,'g'); hold on; 
plot(t,Y-Ysk0,'k');
% plot(t,srr,'r',t,srd_yr+sgyr,'b',t,srd_yr-sgyr,'b'); 
% plot(t,sigsrr+srd_yr,'r:',t,-sigsrr+srd_yr,'r:');
% plot([t(1) t(end)],[srd_yr srd_yr],'r--');
% plot([t(1) t(end)],[srd_yr-sigyr srd_yr-sigyr],'b:' );  
% plot([t(1) t(end)],[srd_yr+sigyr srd_yr+sigyr],'b:' );  hold off
if(~isempty(Yt))
    plot(t,Yt-Ysk0,'b--'); 
    axis('tight'); txalfa='\alpha';
    xlabel(sprintf('Dane v(%d) %s=%.2f',Ldanych,txalfa,alfa));
end
whmin=1.0*stdV; if(whmin<0.1) whmin=0.1; end
whmax=0.8; 
% Etap I. Obliczamy wzmocnienie: 
%     a) ZnaleŸæ Y0 (odpowiedŸ na U0)
%Ypocz=mean(Y(1:nSkok)); %YPocz=mean(Y(1:nSkok-1)); 
%     b) ZnaleŸæ odczinek odpow, gdzie h(t) zmienia siê ma³o: +-5%  
eps=0.05;
okno=round(Tokna/dt); jestUst=0; %*dt
nKpocz=Ldanych-okno+1;
hsr=mean(Y(nKpocz:Ldanych)); jestUst=1;
for(n=1:0) %Ldanych-okno+1:-1:nSkok)
    hsr=mean(Y(n:n+okno-1)); 
    if(max(abs(Y(n:n+okno-1)-hsr))/hsr<eps) if(jestUst==0) nKpocz=n; jestUst=1; end
    else jestUst=0; 
    end
end
Kos=(hsr-Ypocz)/DUskok; % Mamy ju¿ K
Y0lin=Ypocz-Kob*Unom; 
hmin=Kos*DUskok*whmin; hmax=Kos*DUskok*whmax; % zakres wartoœci wyjœcia do identyfikacji

% Etap II. przekszta³camy odpowiedz: % Filtracja: 
%clear y th yh yh1;
y=[]; th=[]; yh=[]; yh1=[];
if(stdV>0.01*Kos)
    Tf=1; oknoF=round(Tf/dt); 
    ncF=round(oknoF/2); 
    np=nSkok-ncF; nf=np+oknoF; % zaczynamy od 
    nk=nKpocz-nSkok; n1=0;
    for(n=n1:nk) 
        nh=n+nSkok; 
        yh(nh)=mean(Y(np:nf))-Ypocz; 
        nf=nf+1; np=np+1; 
    end
    yh(1:nSkok)=0; yh(nh+1:Ldanych)=hsr-Ypocz;
    %zerujemy zaklocenia do chwili skoku 
else %yh(1:nKpocz)=Y(1:nKpocz)-Ypocz; 
    yh=Y-Ypocz; 
end
% yh - przesuniêta do pktu nSkok-1 i wyg³adzona odpowiedz skokowa
Syh=0; nyh=0; lyh=length(yh);
for(i=lyh:-1:10) %round(Tos/dt)) 
    nyh=nyh+1; Syh=Syh+yh(i); if(abs(Syh/nyh/yh(end)-1)<0.01) Ldyh=lyh-nyh-nSkok; end
end
m=0; npocz=-1; nkonc=Ldanych;
% ....... koniec filtracji ...........
Kmaxh=max(yh)/DUskok; wKh=Kos/Kmaxh;
yh1=yh; 
yh=yh1; %*wKh; 
Lyh=length(yh);
for(n=1:Lyh)
    hY=yh(n);
    if(hY>hmin && hY<hmax) 
        m=m+1; y(m)=log(1-hY/(DUskok*Kos)); th(m)=m*dt;% obserwacje wyjœcia
    else
        if(m==0) npocz=n+1; %=n+1; 
        else nkonc=n-1; %=n-1; 
            break; 
        end,
    end
end
if(npocz<-1) npocz=nSkok; fprintf(1,'\nZa silne zaklocenie !!! npocz=%d',npocz); end
subplot(2,1,1); hold on;
if(1)
    plot(([1:Lyh])*dt,yh1+DYshift,'k:',([1:Lyh])*dt,yh+DYshift,'m:'); 
else 
    plot(([nSkok:nSkok+Lyh-1])*dt,yh1+DYshift,'k:',([nSkok:nSkok+Lyh-1])*dt,yh+DYshift,'m:'); 
end
dhmin=Ypocz-Ysk0; 
axis('tight'); ax=axis; 
plot([tSkok tSkok],ax(3:4),'b:',[npocz npocz]*dt,ax(3:4),'m:',...
     [nkonc nkonc]*dt,ax(3:4),'m:',[nKpocz nKpocz]*dt,ax(3:4),'b--'); 
plot([npocz Ldanych]*dt,[hmin hmin]+dhmin,'m:',[npocz Ldanych]*dt,[hmax hmax]+dhmin,'m:')
hold off; 

% Wyznaczamy sta³¹ czasow¹ i opoŸnienie regresj¹ liniow¹
a=polyfit(th,y,1); % a(1)=-1/To; a(2)=d/To; 
Tos=-1/a(1); %/wKh; % To mno¿enie 1/wKh jest zbêdne i blêdne
dels=a(2)*Tos+npocz*dt-tSkok; 
if(dels<=0) dels=dt; end
ds=round(dels/dt); 
t0=nSkok*dt+dels; n0=nSkok+ds;
% .......... Zapis zmierzonej odpowiedzi skokowej ..................
Ldhpom=5*round(Tos/dt); nhPom=nSkok+Ldhpom; if(nhPom>length(Y)) nhPom=length(Y); end
hPom=(Y(nSkok: end)-Ypocz)/DUskok;
% ============= Identyfikacja met.analizy regresji =================
Yref=Y; Ld=Ldanych;  % sygnal referencyjny do selekcji danych
DtR=Tos/10; nDt=round(DtR/dt); % okres interwencji regulatora
skok=1; 
identMNK;
% [Tos,Kos,Kob,DtR,nDt,nKpocz]=IdentOb(Ldanych,t,dt,U,Y,Ysk0,Kob,tSkok,nSkok,Unom,Ypocz,Yt,stdV,DUskok,DYshift)