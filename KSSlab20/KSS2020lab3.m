% Ekonometria cwicz.1
clear all; 
global Unom FU0 X0 Y0lin rzad ToSum Tob Kob K1 K2 dt invh DtR WspXn Nielin C1 Cf  S1 Sf;
global hPom hEmp alfa; 
global przyr wagi KSSlab PokazId delr txPar; 
KSSlab=1;
przyr=1; przyr=0; wagi=1; PokazId=1; 
% ======= Wybór typu obiektu ============================
invh=1; % invh=0 obiekt nieminimalnofazowy (trudny); =0 minimalnofazowy (³atwiejszy)  
Nielin=1; % Nielin=0 obiekt liniowy; =1 nieliniowy bez interakcji poziomów; =2 - z interakcj¹ poziomów; =2;
% ................. Podstawienia parametrów symulacji .....................
dt=0.001; Tsym=12; if(Nielin>1) Tsym=21; end
Ldanych=round(Tsym/dt); Lbins=30; 
Nsym=Ldanych; %round(Tsym/dt);
% teoretyczne param.rozkladow generowanych
srtv=0; sigw=0.01; srtw=0.0; 
Sigv=0.1; %0.03; %0.1;%.05; %.05; %1; %0.01; %Sigv=0.05; %Sigv=0.1; 
% Dane dla zaklócen
Tau=0.1; %dt*10; %100;  
Tokna=2; okno=round(Tokna/dt); jestUst=0; %*dt
nSkok=round(Ldanych*1/3); tSkok=nSkok*dt; 
alfw=.98; 
alfa=exp(-dt/Tau); Kw=sqrt(1-alfw^2); 
% ......... Podstawienia startowe symulatora zak³óceñ ARMA
yARv1=0; yARv2=0; yARw=0; 
brn=0; dU=0; Up=1;
% =================== Dane dla procesu ==================
rzad=5; ToSum=1; Ts=ToSum; Tob=Ts/rzad; Kob=2; Tob2=Tob;
% Obiekt nieliniowy 1 (Nielin=1) bez interakcji poziomów: 
% .......................................................
% Na wejœciu mamy przep³yw Fwe, zasilaj¹cy ci¹g zbiorników, zbocznikowany dodatkowym zbiornikiem. 
% Wyjœciem jest przep³yw Fwy, bêd¹cy sum¹ przep³ywów Fwy1 z ci¹gu g³ównego i Fwy2 z bocznika. 
% Celem regulacji jest reakcja na zmianê przep³ywu Fwe, tak aby jak najszybciej osi¹gn¹æ Fwy=Fwe,
%   a tak¿e jak najszybciej uzyskaæ ustalony poziom w ci¹gu zbiorników (poziomy s¹ wektorem stanu X). 
% Osi¹ga siê to zmieniaj¹c wartoœæ zadan¹ wyjœcia Yref dla Fwy, 
%   do wartoœci zmierzonego przep³ywu Fwe. 
% Sterowaniem U jest ustawienie zasuwy na wejœciu, rozdzielaj¹cej strumieñ Fwe na
%   uFwe1 do ci¹gu g³ównego i uFwe2 do zbiornika bocznikuj¹cego, przy sta³ym ciœnieniu P0 przed zasuw¹.  
% OdpowiedŸ skokow¹ Y=Fwy na zmianê U pokazuje Figura 1, a odpowiedzi skokowe
%   poziomów Poz(1:rzad) oraz sk³adowych Fw1 i Fw2 na zmianê F1 i F2 wynikaj¹c¹ ze skoku U
%   przedstawia Figura 13.
% ...................................................................
%  Fwe 
% ==X===>|~~~~~~|
%   |    | Zb.1 |_
%   |    \_______=>Fwy1
%   }            |~~~~~~|
%   | Fbocz      | Zb.2 |
%   | dla invh=1 \_______=>Fwy2
% ..|................................... =====>Fwy_rzad-1
%   | dla invh=0 Fbocz=0                     |~~~~~~~|
%   |                                        | Zb.rz |
%   |                                        \________=>Fwyrzad
%   |==================>=======>Fbocz                  |
%                          |~~~~~~~~~~~|               | 
%                          | Zb.boczn  |               |
%                          \____________=>Fyb==========+===> Fwy

% Obiekt nieliniowy 2 (Nielin=2) z interakcjami poziomów: 
% ==X====> Fwe
%   |    |~~~~~~|
%   |    | Zb.1 ||~~~~~~|
%   v    |      || Zb.2 ||~~~~~~|
%   }    |      ||      || Zb.3 |
%   |    |      ||      ||      ||~~~~~~| 
% ..|...................||...................|............    
%   |    |      ||      ||      || Zb.4      ||~~~~~~~~~|
%   |    |      ||      ||      ||           || Zb.rzad | +
%   |    |______________________________________________==>X===>Fwy
%   |==================>=======>Fbocz                      |+
%      dla invh=0 Fbocz=0       |~~~~~~~~~~~|              ^ 
%                               | Zb.boczn  |              |
%                               \____________=>Fyb========>| 
% =============================================================================== 
Umin=0; Umax=2; % ograniczenia sterowania
Unom=1; Up=Unom;
if(invh) 
    % ......... Obiekt nieminimalnofazowy - dobór parametrów bocznika ....
    aa=0.5; aa=0.1; aa=0.2;
    aw=(1-aa+aa*(rzad-1))
    wK2=aw/(rzad-1); %aw>1 i aw<(rzad-1)
    K1=Kob/(1-wK2); 
    Ko=K1^(1/(rzad-1)); K2=wK2*K1;
    wKT=wK2*Tob2/Tob; T2=Tob2;
else K1=Kob; K2=0; wKT=1; % obiekt minimalnofazowy
end
% Stan poczatkowy i parametry: 
clear Y;
Fnom=10; Y(1:rzad)=Fnom; 
if(Nielin) 
    [F1nom,F2nom,Fsum,X0]=parProc(Unom,Kob,wKT,Ts,Tob,Tob2,Nielin,invh);
    X=X0; Y0lin=Fnom-FU0; 
else
    % Obiekt liniowy - szereg inercji I.go rzêdu 
    %    dla invh>0 dodatkowo bocznik, tj. jedna inecja po³¹czona równolegle
    WspXn=-1; rz1=rzad-invh; u0=0.2; u2=(1-u0)/2; 
    X(1)=K1*Up; for(i=2:rz1) X(i)=X(i-1); end 
    if(invh) 
        X(rzad)=K2*Up; Fsum=X(rz1)-X(rzad); 
    else Fsum=X(rzad); 
    end
    X0=X; Y0lin=Fnom-Kob*Unom; FU0=0; 
end
% ........................................
srtv=Fnom-Fsum; 
sigv=Sigv*Kob*Unom; 
Kv=sigv*sqrt(((1-alfa^2)^3)/(1+alfa^2)); 
% =============== Petla RT =====================
clear y U Yt v;
ntr=1; v(1:ntr)=srtv; t(1:ntr)=[0:ntr-1]*dt; Yt(1:ntr)=Fnom; Poz(:,1:ntr)=X; 
Ystart=0;
for(n=1:Ldanych)
    t(n+ntr)=(n+ntr)*dt;
    U(n)=Up+dU;
    % ------ Zaklocenie wyjsciowe ---------------
    zv=Kv*randn; % rozkl.normalny jest juz standaryzowany
    yARv1=yARv1*alfa+zv; yARv=yARv1;  
    yARv2=yARv2*alfa+yARv1; 
    yARv=yARv2; 
    v(n+ntr)=yARv+srtv;
    % ------ Zaklocenie wejsciowe ---------------
    zw=randn; % rozkl.normalny jest juz standaryzowany
    zw=zw*sigw; 
    yARw=yARw*alfw+Kw*zw;
    w(n)=yARw+srtw;
    % Symulacja obiektu: X zmienne stanu Y wyjœcie
    [X,y,uu,uF1(n),uF2(n),yFwy1(n+ntr),yFwy2(n+ntr)]=Proces(X,U(n),0); 
    Yt(n+ntr)=y+srtv; %FU0;
    Poz(:,n+ntr)=X; 
    Y(n+ntr)=y+v(n+ntr); % wyjscie 
    % ..... Sterowanie skokowe do identyfikacji czynnej .........
    Up=U(n); 
    if(n==nSkok-1) 
        Ypocz=mean(Y(1:nSkok-1)); stdV=std(Y(1:nSkok-1)); 
        Ystart=Ypocz-Kob*Unom; 
        if(stdV<0.01) dU=0.05; 
        else dU=5*stdV; if(Unom+dU>Umax) dU=Umax-Unom; end,  %dU=(Umax-Unom); 
        end
        DUskok=dU; 
    else dU=0; 
    end
    % ...........................................................
end
% Start symulacji odpowiedzi skokowej ........................
DYshift=1+1.5*DUskok; % przesuniêcie poprawiaj¹ce widocznoœæ
Ysk0=Ypocz-DYshift; 
Yt(1:ntr)=Yt(1+ntr); % Yt=Yt-Ysk0; %Ypocz+DYshift; %Yt=Yt-Ystart;
U(n+1:n+ntr)=U(n); Upocz= U(1); 
vh=v;  
figure(1); 
subplot(2,1,2); 
plot(t(1:end-ntr),w,'r:',t,v,'k'); %hold on; 
axis('tight');  
txalfa='\alpha';
xlabel(sprintf('Dane w(%d) r; v(t) k: %s=%.2f',Ldanych,txalfa,alfw));
% ............. Wykresy poziomów ....................
figure(13); kol='kbmrgckbmrgckbmrgc'; if(invh) kol(rzad)='k'; end
subplot(2,1,1);
for(r=1:rzad) plot(t(1:end-ntr),Poz(r,ntr+1:end),kol(r)); hold on; end
axis('tight'); ax=axis; Dax=ax(4)-ax(3); axis([ax(1:2) ax(3)-0.05*Dax ax(4)+0.05*Dax]);
xlabel('Odpowiedz skokowa poziomów w zbiornikach');  hold off,
subplot(2,1,2);
plot(t(1:end-ntr),uF1(1:end),'k-.',t(1:end-ntr),uF2(1:end),'r-.',t(1:end-ntr),uF1+uF2,'m:',...
     t(1:end-ntr),yFwy1(ntr+1:end),'k',t(1:end-ntr),yFwy2(ntr+1:end),'r',t(1:end-ntr),yFwy1(ntr+1:end)+yFwy2(ntr+1:end),'m');
axis('tight'); ax=axis; Dax=ax(4)-ax(3); axis([ax(1:2) ax(3)-0.05*Dax ax(4)+0.05*Dax]);
xlabel('Wejscia kanalow F1 k-.; F2 r-.; F=F1+F2 b.-, i wyjscia Fw1 k; Fw2 r; Fw=Fw1+Fw2 g'); 
% ............. Identyfikacja obiektu G(s)=K*exp(-d)/(1+Ts) .............
clear w; w=0;
 [hPom,Tos,Kos,dels,ds,Tor,Kor,delR,dr,Kob,DtR,nDt,nKpocz]=IdentOb(Ldanych,t,dt,U,Y,Ysk0,tSkok,nSkok,Unom,Upocz,Ypocz,Yt,stdV,DUskok,DYshift,ntr,rzad,Kob,Ts,Tokna);       
% ======= Synteza regulat. predykcyjnego DMC ======
Nu=5; wR=1;
nSkok=round(Ldanych*1/6); tSkok=nSkok*dt; 
% ================= Regulacja ===================== 
Ld=1500; %0; 
txDel='\Delta'; txthet='\theta'; 
% .......... Generujemy profil wart.zadane Yref .......................
clear Yref; 
TsRef=1; tS=round(TsRef/dt); tS1=15*tS; 
Ld=30*tS; %round(Ld*2/3); 
if(Nielin && invh) d1Yref=0.2*(Y(1)-Y0lin); d2Yref=-0.2*(Y(1)-Y0lin);
else d1Yref=0.6*(Y(1)-Y0lin); d2Yref=-0.7*(Y(1)-Y0lin);
end
Yref(1:tS)=Y(1);Yref(tS+1:tS1)=Yref(tS)+d1Yref; Yref(tS1+1:Ld+ntr)=Yref(tS)+d2Yref; 
% .............. Ograniczenia dla regulatora ..........................
DU_tech=1*dt; %0.5*dt; % ograniczenie przyrostu sterowania (techniczne)
dUtech=DU_tech*1000; % likwid. ograniczenia techniczn.
DUreg=DU_tech/2*nDt; % ograniczenie dla regulat. na okres regulacji
% ............... Symulacja regulacji .......................
for(nrSym=1:1)
    RegulPID;
    RegulPred; if(nrSym==2) break; end
% =================== Identyfikacja bierna ==============================
    DtReg=DtR; nDtR=nDt;
    %KSSident;
end


