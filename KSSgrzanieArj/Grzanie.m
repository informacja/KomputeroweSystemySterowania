%Grzanie
charKotla=1;
if(charKotla==1)
    clear all
    global dt wExp KqSciIz KqScian KqScIzSum KTq Tauw Mpowcw MSccw Kkal CF Cq Kkal_Cq cicho ro0 cw;
    global juzByleta Tob Kob rzads rzad rps KTTs KTq0;
    global wAqTqT Aq0 AqTg AqT KTsT KTsTzw Taus Mpow  Fvwz Fwconst;
    global ndrg ndrg;
    global przyr wagi nfig1 invh;
    global Unom X0;
    global hPom hEmp nTwg nTwg1 lX KSSlab PokazId;
    KSSlab=0; PokazId=1;
    przyr=1; przyr=0; wagi=1; nfig1=1; invh=-1;
    nTwg=9; nTwg1=nTwg+1; lX=4+nTwg;
    charKotla=1;
    TzewSr=-10; Tnom=20; T=Tnom; % nominalna temperat. w pomieszcz.
    
    % Dobowy profil temperat. od 6.00 do 24.00 co 6 h
    dt=1/3600; dt=1.e-3;
    Tzw=[-4  2  0 -3]; lt=length(Tzw);
    tg= [ 6 12 16 24]; a0=mean(Tzw); ntsin=lt/2;
    for(k=1:lt)
        if(k<=ntsin) Fi(1:lt,k)=sin(2*pi*k/24*tg); %Fi(1:lt,2)=sin(2*pi/12*tg); Fi(1:lt,3)=sin(2*pi/8*tg);
        else    Fi(1:lt,k)=cos(2*pi*(k-ntsin)/24*tg); %Fi(1:lt,5)=cos(2*pi/12*tg); Fi(1:lt,6)=cos(2*pi/8*tg);
        end
    end
    Ar=inv(Fi'*Fi)*Fi'*(Tzw-a0)'; %A=fft(Tzw-a0);
    Tr=Fi*Ar+a0;
    %A=fft(Tz-a0); af0=A(1); as(1)=real(A(5)); ac(1)=imag(A(5)); as(2)=real(A(6)); ac(2)=imag(A(6));
    Tsym=48;
    t=[6:dt:Tsym+6]; Lt=length(t);
    Tz(1:Lt)=0; TzewSr;
    for(k=1:lt)
        if(k<=ntsin) Tz=Tz+Ar(k)*sin(2*pi*k/24*t); else Tz=Tz+Ar(k)*cos(2*pi*(k-ntsin)/24*t); end, %sin(2*pi/12*t)+Ar(2)*sin(2*pi/12*t)+Ar(3)*sin(2*pi/8*t)+...
        %Ar(4)*cos(2*pi/24*t)+Ar(5)*cos(2*pi/12*t)+Ar(6)*cos(2*pi/8*t);
    end
    Tzewn=Tz-mean(Tz)+TzewSr; clear Tz;
    A0=Tzewn(1)-Tzw(1);
    figure(101);
    plot([tg 24+tg],[Tzw Tzw]+A0,'k.',[tg 24+tg],[Tr;Tr]+A0,'ro',t,Tzewn,'b'); axis('tight'); grid;
    xlabel(sprintf('Profil dobowy temperatury zewnêtrznej: TzewSr=%.2f ^oC',TzewSr));
    
    % ===================== Szacunki zapotrzebowania na cieplo ===============
    % https://muratordom.pl/instalacje/ogrzewanie-domu/dobieramy-moc-grzewcza-pieca-do-co-poradnik-aa-EoA1-uSAM-SGer.html
    %     od 120 do 200 W/m? – dla domów bez izolacji cieplnej, wybudowanych przed rokiem 1982;
    %     od 90 do 120 W/m? – dla domów z lat 80. i 90. ubieg³ego wieku;
    %     od 60 do 90 W/m? – dla domów wznoszonych od koñca lat 90. ubieg³ego wieku, dobrze zaizolowanych, z nowoczesnymi oknami;
    %     50 W/m? -  dla domów budowanych zgodnie z wymaganiami obowi¹zuj¹cymi od 2017 r.
    %     W nowym domu o powierzchni 150 m2 zapotrzebowanie na ciep³o wynosi mniej wiêcej 10 000 W, czyli 10 kW.
    %     W ofercie producentów s¹ kot³y maj¹ce moc nie mniejsz¹ ni¿ 12-15 kW.
    % Wzap=50*120;
    WzapMax=10; Wzap=WzapMax/1.4;
    qzap=Wzap*3600; % Wzap [kW] qzap [kW/h]
    Lc=0.25; Liz=0.2; % Lc gruboœæ sciany ceglanej; Liz grub.izolacji;
    Spom=121; Hpom=2.8; z=sqrt(Spom); s=4*z; % s=44m;
    L=Lc+Liz;
    % ......... Stale fizyczne ..........
    cw=4.19; % cieplo wlaœciwe [kJ/kg]
    wExp=0.21e-3; hrg=4; % wspolcz. rozszerzalnoœci wody i dlugoœæ ci¹gu
    ro0=1000/(1+wExp*16); % gest. w 20. stC
    % ..................................................
    Fwconst=1; 
    Twg=70; Twz=50; % Szacunkowe temp. wody goracej i zimnej
    Twg=45; Twz=25; 
    %Twg=55; Twz=50; % Szacunkowe temp. wody goracej i zimnej
    roz=ro0/(1+wExp*(Twz-20)); rog=ro0/(1+wExp*(Twg-20));
    Fw=qzap/(cw*(Twg-Twz)); % sredni przepl. masowy Fw=601 kg/h=0.167 kg/s=0.170 l/s; % kg/h  
    Fvwz=qzap/(cw*(rog*Twg-Twz*roz)); % Faktyczny (sta³y) przeplyw wody zimnej (praca pompki pieca)
    % ............................
    %Fvwz=Fwz/roz; Fwg=Fvwz*rog; 
    qFw=cw*Fvwz*(rog*Twg-Twz*roz);
    %qFw-qzap,
    % Cisnienie ciagu grawit DPg=wExp*hrg*rog*(Twg-Twz)=Fw*mi_sr*RrK;
    DPg=wExp*hrg*rog*(Twg-Twz); % (16.6 Pa) kg*m/s^2/m^2;
    % DPg=CF*RrK*rog*(Twg-Twz)=Fw*mi_sr*RrK;
    juzByleta=0;
    [mi_sr miw etaw]=lepkw(Twg,Twz,T);
    % DPg=RrK*eta*Fw/rog
    % DPg=RrK*mi_sr*Fw
    RrK=DPg/(mi_sr*Fw); %*100; F=DP/(mi_sr*RrK); F=(DP/DPg)*mi_sr*Fw
    % CF - wspolczynnik dla Fw=CF*rog/misr*(Twg-Twz)
    CF=wExp*hrg/RrK;
    Cq=CF*cw;
    % Straty ciepla qstr=Kqz*(Ts-Tzewn)=KqScian*(T-Ts);
    % .......... Obliczenie parametrów cieplnych budynku ...................
    %       qzap-Fp*cwp*(T-TzewSr)-(T-Ts)*KqScian=0 - bilans w pomiesz.
    %       (T-Ts)*KqScian-(Ts-TzewSr)*KqSciIz=0    - bilans w scianie
    %       Ts*(KqScian+KqSciIz)=T*KqScian+TzewSr*KqSciIz
    % Przyjmujemy KqScian=wK*KqSciIz; np. wk=9;
    %       (T-Ts)*wK*KqSciIz-(Ts-TzewSr)*KqSciIz=0
    %       T*wK=Ts*(1+wK)-TzewSr
    %       qzap-Fp*cwp*(T-TzewSr)-(T-Ts)*KqScian=0
    T=Tnom; wK=10; %wK przelicznik przewodnoœci ceg³a/styropian =0.72/0.045=16;;
    % Oszacowanie Ts - sredniej temp. œcian 
    Ts=(T*wK+TzewSr)/(1+wK); %Ts*(KqScian+KqSciIz)=T*KqScian+TzewSr*KqSciIz
    % ...... Parametry cieplne budynku ..........
    % Straty przez przegrody daj¹ 70% ca³ego zapotrzebowania,
    % a 30% to wentylacja Fp:
    cwp=1.; %2.02; % cieplo w³aœciwe powietrza przy stalym ciœn. i T=20 stC [kJ/kg]
    ropow=1.12; Hpom=2.8; Mpow=(Hpom*Spom)*ropow;
    % ....................................
    Fp=0.3*qzap/(cwp*(T-TzewSr)); % Oszacowanie przeplywu wentylacji (tak, aby zabiera³a 30% ciep³a) 
    % Oszacowanie sta³ej przewodnoœci scian
    KqScian=(qzap-Fp*cwp*(T-TzewSr))/(T-Ts); KqSciIz=KqScian/wK;
    KqScIzSum=(KqSciIz+KqScian);
    % ...... Wzmocnienia równan dynamiki
    %KqScian - sta³a cieplna scian; KqScIzSum - œrednia sta³a œcian KTq -
    KTq=1./KqScian; KTsT=KqScian/KqScIzSum; KTsTzw=KqSciIz/KqScIzSum;
    % ...... Parametry geometryczne budynku
    Spom=120; %[50 25 20 15 10]; % powierzcnie
    Ssum=sum(Spom); uSpom=Spom/Ssum;% powierzchnie pomieszczeñ razem 120
    % Teraz parametry dynamiki:
    % Tauw=Mpow*cw/KqScian/3600; Mpowcw=Vpow*ropow*cw; ropow=29*1.e5/(8310*293)=1.12 kg/m^3
    Tauw=Mpow*cwp/KqScian; %/3600 %Tauw=[0.088 0.0352 0.0264 0.0264] h
    Mpowcw=Mpow*cwp;
    Taus=5; % szacunkowa sta³a czasowa scian
    MSccw=Taus*KqScIzSum; Ms=MSccw/0.80; % 0.8 szac.cw scian;  Ms=40*2.8*2*0.45*1.8e3=181000, a Ms =189000
    Taufp=Mpow/Fp;
    % ----------- Rownania dynamiki ogrzewania - parametry
    %qzap=Cq*ro^2(Twg-T)^2*(1-exp(- (Kkal/Cq)/(ro*(Twg-Twz)) ))^2
    %sqrt(qzap/Cq)/(ro*(Twg-Tnom))=(1-exp(-(Kkal/(CF*ro^2))/(Twg-Twz)));
    %(Kkal/(CF*ro^2))/(Twg-Twz)=-log(1-sqrt(qzap/Cq)/(ro*(Twg-Tnom)));

    % ========== Model stanu ustalonego kaloryfera: we-wy kaloryfera: 
    % Wzór na róznicê temperatur wody: 
    % Twg-Twz=(Twg-T)*(1-exp(-Kkal/(Fw*cw))); Fw - œredni przeplyw masowy
    % Twz-T=(Twg-T)*exp(-Kkal/(Fw*cw))); 
    % (Twg-T)=(Twz-T)*exp(Kkal/(Fw*cw))); 
    % St¹d: (1-(Twg-Twz)/(Twg-T))=exp(-Kkal/(Fw*cw))
    % Obliczamy wymagan¹ sta³¹ Kkal: 
    % Kkal/(Fw*cw)=-log(1-(Twg-Twz)/(Twg-T));  
    Kkal=-Fw*cw*log(1-(Twg-Twz)/(Twg-T)); 
    Kkal_Cq=(Kkal/Cq);
    DTn=(Twg-T)*(1-exp(-Kkal_Cq*CF/(Fw))); wx=[DTn (Twg-Twz)/(Twg-T)]
    Fn=CF*rog/mi_sr*DTn; %*(1-exp(-Kkal/(Fw*cw)));
    %Kkal=-log(1-sqrt(qzap/Cq)/(rog*(Twg-Tnom)))*(Twg-Twz)*CF*rog^2;
    
    qKal=Fn*cw*DTn; %wq=qKal/(Cq*rog/mi_sr*DTn^2)-1,; %wq=qKal/qzap-1,
    qkal=qzap; Twgp=Twg; Twzp=Twz; Fwp=Fw; cicho=-1;
    %Kociol;
    % gzap>0 charkterystyka qwy(Twg(1:lwg));
    % gzap<0 obliczenie qwy(Twg);
    if(Fwconst)
        figure(202);
        Ldq=100; dlaFwy=0;
        if(dlaFwy)
            dTg=(90-Twz)/(Ldq-1); Fw0=Fw;
            Twgd=[Twz:dTg:90]; 
            for(i=1:Ldq)
                %[Twgn(i),qKaln(i),brakTwg(i),rog,roz]=KociolFwz(T,Twg,Twz,qZap(i),Fwz);
                [y,Twgy(i),Twg1,qwy(i),Fwe(i),Fwy(i),rog]=modKotla(T,Twgd(i),Twz,qzap,Fwz,roz);
            end
            subplot(2,2,1); plot(Twgd,Twgy,'k'); axis('tight');
            subplot(2,2,3);
            plot(Twgd,Fwe,'k',Twgd,Fwy,'r',...
                [Twgd(1) Twgd(end)],[Fwz Fwz],'m',Twg,Fw0,'r*'); axis('tight');
            subplot(1,2,2); plot(Twgd,qwy,'k',[Twgd(1) Twgd(end)],[qzap qzap],'r'); axis('tight');
            
        else
            dq=0.8; dqZ=dq/(Ldq-1); Twg0=Twg; Twz0=Twz;            
            qZap=[1-dq/2:dqZ:1+dq/2]*qzap;
            Tx=[14:3:26]; lTx=length(Tx);  %Tx=[14:2:26]; lTx=length(Tx);           
            yTz=Tx; 
            for(k=1:lTx)
                for(i=1:Ldq)
                    [Twgn(i),Twzn(i),qKaln(i),Fwn(i),Fwy(i),brakTwg(i),rog,roz]=KociolFwz(Tx(k),Twg,Twz,qZap(i));
                    Twg=Twgn(i); Twz=Twzn(i);
                end
                qxZap=qZap/3600; qxzap=qzap/3600; qxKaln=qKaln/3600;
                yFw=Fwn; 
                subplot(2,2,1); plot(qxZap,qxKaln,'k',qxZap,qxZap,'r:');
                xlabel('qKal,qZap [kWh]'); axis('tight'); hold on;
                %subplot(2,2,3); plot(qxZap,brakTwg,'k'); xlabel('brak'); axis('tight');
                subplot(2,2,4); plot(qxZap,Twgn-Twzn,'k');
                xlabel('Twg-Twz'); axis('tight');  hold on;
                subplot(2,2,2);
                plot(qxZap,Twgn,'r',qxZap,Twzn,'b',qxzap,Twg0,'r*',...
                    qxzap,Twz0,'b*'); xlabel('Twg,Twz'); axis('tight');  hold on;
                subplot(2,2,3); plot(qxZap,Fwn,'r',qxZap,Fwy,'b--');
                hold on; xlabel(sprintf('Fwn,Fwy [kg/h]; Fvwz=%.3g [dcm^3/s]',Fvwz/3.6)); axis('tight');
                % qKal=Aq0+AqTq*Twg+AqT*T; 
                Fiq=[]; Fiq(1:Ldq,1)=1; Fiq(1:Ldq,2)=Twgn'; 
                G=Fiq'*Fiq; AqTgx(:,k)=inv(G)*Fiq'*qxZap'; 
                xTwgn(:,k)=Twgn'; 
            end
            for(k=1:4) subplot(2,2,k); hold off; end
            Txgp=min(Twgn); Txgk=max(Twgn); Txzp=min(Twzn); Txzk=max(Twzn); 
            Txzp=14; Txzk=28;             
            Fiq=[]; Fiq(1:lTx,1)=1; Fiq(1:lTx,2)=Tx'; 
            G=Fiq'*Fiq; AqT2=inv(G)*Fiq'*AqTgx(1,:)'; 
            Aq0=AqT2(1); AqTg=mean(AqTgx(2,:)); AqT=AqT2(2); wAqTqT=0;
            figure(201); kol='mrbkgc'; kk=0; tx='qKal=f(Tg,T); '; 
            %Tx=[14:3:26]; lTx=length(Tx);  
            for(k=1:lTx)
                qK=[]; 
                qK=Aq0+AqTg*xTwgn(:,k)+AqT*Tx(k); 
                kk=kk+1; 
                plot(xTwgn(:,k),qK,kol(kk),xTwgn(:,k),qxZap,'k--'); hold on; %[kol(kk) '--']); hold on; 
                tx=[tx sprintf('T_%c=%.1f; ',kol(kk),Tx(k))]; 
            end
            hold off; axis('tight'); xlabel(tx);
            % qKal=Aq0+AqTq*Twg+AqT*T; 
        end
    else
        [Fwn,Twzn,qKal,Twgn]=Kociol(T,Twg,Twz,qzap);
        wx=[Fwn Twz Twg;Fwp Twzp Twgp], 
        Txgp=40; Txgk=80; Txzp=12; Txzk=34; 
    end
    %[qKal qzap],
    %return;
end % sztuczny koniec
if(charKotla>=0 && Fwconst==0)
    % ========== obliczymy charakterystyké kotla Tg=f(q,T) =====
    cicho=1;
    Txg=[Twgp:1:Txgk]; lTg=length(Txg); lTx=40;
    Txx=rand(lTx,1)*(30-10)+10; % To musi byc zrandomizowane
    Tx=sort(Txx); clear xqK xTx yTg Txx G Fi np nk;
    Tx=[Txzp:2:Txzk]; lTx=length(Tx); 
    Twzy=Twz;
    n=0;  
    for(i=1:lTx)
        np(i)=n+1;
        for(k=1:lTg)
            if(Fwconst)
                [Twgn,Txzn,xq,Fwx,Fwy,brak]=KociolFwz(Tx(i),Txg(k),Twzy,qzap);
            else [Fwx,Txzn,xq,Twgn,brak]=Kociol(Tx(i),Txg(k),Twzy,-qzap);
            end
            if(brak==0) 
                n=n+1; 
                xTx(n)=Tx(i); xTg(n)=Txg(k); 
                yqK(n)=xq/3600; Twzy=Txzn; 
                yFw(n)=Fwx; yTz(n)=Txzn; 
            end
        end
        nk(i)=n;
    end
    x1=yqK; y=xTg; % y=yqK; x1=xTg;
    Fi(1:n,1)=1; Fi(1:n,2)=x1(1:n)'; Fi(1:n,3)=xTx(1:n)'; 
    % Tg=A1+A2*qK+A3*T;
    if(Fwconst==0)
        cztery=4;
        if(cztery) Fi(1:n,4)=(xTx(1:n)').*(x1(1:n)')/10; end
        [V,S,P]=svd(Fi'*Fi);
        U=Fi;%*P;
        G=(U(1:n,:)'*U(1:n,:)); ATwg=inv(G)*U(1:n,:)'*y(1:n)'; 
        if(cztery) ATwg(4,1)=ATwg(4,1)/10; end
    else
        cztery=0; 
        G=Fi'*Fi; ATwg=inv(G)*Fi'*y(1:n)'; 
        % Twg=ATwg(1)+ATwg(2)*qKal+ATwg(3)*T;
    end
    % --------------- Przeliczamy model na q=f(Tg,T) -------------
    Aq0=-ATwg(1)/ATwg(2); AqTg=1/ATwg(2); AqT=-ATwg(3)/ATwg(2); 
    if(cztery) wAqTqT=ATwg(4)/ATwg(2); else wAqTqT=0; end
    % ........ Obliczamy model .....
    Txs=Tx; lTxs=lTx; Tx=[Txzp:2:Txzk];  lTx=length(Tx); 
    kol='kbrmgc'; txp='Zalezn. Twg(T,q):'; txpl='T=';
    j=1;
    figure(202)
    for(i=1:lTx)
        %while(j<lTxs && Txs(j)<Tx(i)) j=j+1; end
        xQ=[]; yTwg=[];
        xqmin=min(yqK(np(i):nk(i))); xqmax=max(yqK(np(i):nk(i)));
        dq=(xqmax-xqmin)/(lTg-1); xQ=[xqmin:dq:xqmax];
        clear Fi; 
        for(k=1:lTg)
            Fi(1:lTg,1)=1; Fi(1:lTg,2)=xQ'; Fi(1:lTg,3)=Tx(i); 
            if(cztery) Fi(1:lTg,4)=Tx(i)*(xQ'); end
        end
        yTg=Fi*ATwg;
        % Yg=f(q,T)
        YqK=(Aq0+AqTg*Txg+AqT*Tx(i))./(1+wAqTqT*Tx(i)); 
        %Aq0=-ATwg(1)/ATwg(2); AqTg=1/ATwg(2); AqT=-ATwg(3)/ATwg(2); 
        % Twg=ATwg(1)+ATwg(2)*qKal+ATwg(3)*T;
        % qKal=Twg/ATwg(2)-ATwg(1)/ATwg(2)-ATwg(3)/ATwg(2)*T; 
        % qKal=Aq0+AqTq*Twg+AqT*T; 
        if(i<7) 
             Kol=[kol(i) ':']; Koly=kol(i); 
        else Kol=[kol(i-6) ':']; Koly=[kol(i-6) '--']; 
        end
        txpl=[txpl sprintf('%.0f %s; ',Tx(i),Kol)];
        subplot(2,2,1)
        plot(yqK(np(i):nk(i)),xTg(np(i):nk(i)),Koly,xQ,yTg,Kol); hold on;
        hold on;
        subplot(2,2,3);  
        plot(xTg(np(i):nk(i)),yqK(np(i):nk(i)),Koly,Txg,YqK,Kol); hold on; 
        subplot(2,2,2) 
        plot(yqK(np(i):nk(i)),yFw(np(i):nk(i))/3600,Koly); hold on;
        hold on;
        subplot(2,2,4);  
        plot(xTg(np(i):nk(i)),yTz(np(i):nk(i)),Koly); hold on;         
    end
    subplot(2,2,1); axis('tight'); xlabel([txp txpl '    q [kWh]']); ylabel('Twg ^oC'); hold off; 
    txp='Zalezn. q=F(Tg,T)'; 
    subplot(2,2,3); hold off; axis('tight'); xlabel([txp txpl '    Twg ^oC']); ylabel('q [kWh]'); 
    subplot(2,2,2); axis('tight'); ylabel('Fw [kg/s]')
    subplot(2,2,4); axis('tight'); ylabel('Twz [^oC]')
%    return;
end
% ---------- Przeliczenie parametrow dynamiki ---------------
KTTgp=(KTq*AqTg)*3600; KTT=(KTq*AqT)*3600; KT0=KTq*Aq0*3600;
Tob=Tauw/(1-KTT); KTTs=1/(1-KTT); Kob=KTTgp/(1-KTT); KTq0=KT0/(1-KTT); 
% ...........................................................
% .... Stan ustalony qK i Ts dla zad. T i Tzewn .....
% T/Tauw=Ts*(KqScian/Mpowcw)+qK*(1/Mpowcw);
% Ts/Taus=(T*KqScian/MSccw+Tzewn*KqSciIz/MSccw);
% Dla pomieszcz.:
% 0=-(Fp/Mpow)*(T-Tzewn)-(T-Ts)*KqScian/Mpowcw+qK/Mpowcw;
% 0=-T*(Fp/Mpow+KqScian/Mpowcw)+Tzewn*(Fp/Mpow)+Ts*KqScian/Mpowcw+qK/Mpowcw;
% qK=Aq0*3600/(1+wAqTqT*T)+AqTg*3600/(1+wAqTqT*T)*Twg+AqT*3600/(1+wAqTqT*T)*T;
% Parametry dynamiki dla u=Twg:
% dT_dt=(-T*(Fp*cwp+KqScian-AqT*3600/(1+wAqTqT*T))+Tzewn*(Fp*cwp)+Ts*KqScian+qKr)/Mpowcw;
% Taupom=Mpowcw/(Fp*cwp+KqScian-AqT*3600/(1+wAqTqT*T));
% dT_dt=-T/Taupom+(Tzewn*(Fp*cwp)+Ts*KqScian+qK0+AqTg*3600/(1+wAqTqT*T)*Twg)/Mpowcw;
% Kob=AqTg*3600/(1+wAqTqT*T)/Mpowcw*Taupom;

% T=(Tzewn*(Fp/Mpow)+Ts*KqScian/Mpowcw+qK/Mpowc)/(Fp/Mpow+KqScian/Mpowcw)
% lub:
% qK/Mpowcw=T*(Fp/Mpow+KqScian/Mpowcw)-Tzewn*(Fp/Mpow)-Ts*KqScian/Mpowcw;
% qK=T*(Fp*cwp+KqScian)-Tzewn*(Fp*cwp)-Ts*KqScian;
% Dla sciany: 
%(T-Ts)*KqScian=(Ts-Tzewn(1))*KqSciIz
%T*KqScian+Tzewn(1)*KqSciIz=Ts*(KqScian+KqSciIz)
% Parametry dynamiki T->Twg
Taupom=Mpowcw/(Fp*cwp+KqScian-AqT*3600/(1+wAqTqT*T));
Kob=AqTg*3600/(1+wAqTqT*T)/Mpowcw*Taupom;
% Stan ustalony
Ts=(T*KqScian+Tzewn(1)*KqSciIz)/KqScIzSum; 
%Tss=Taus*(T*KqScian+Tzewn(1)*KqSciIz)/MSccw;
Twgx=Twg0; Twzx=Twz0;
qK=T*(Fp*cwp+KqScian)-Tzewn(1)*(Fp*cwp)-Ts*KqScian;
qKs=Mpowcw*Tnom./Tauw-Ts*KqScian; %qK=4.e4; %qK=20.e4;
cicho=1;
if(Fwconst==0)
    [Fwn,Twzn,qKaln,Twgn,brak,rog,roz]=Kociol(T,Twg,Twz,qK);
    %tx='[Fwn,Twzn,qKaln,Twgn]./[Fw,Twz,qK,Twg]', [Fwn,Twzn,qKaln,Twgn]./[Fw,Twz,qK,Twg],
else
    [Twgn,Twzn,qKaln,Fwn,Fwyx,brak,rog,roz]=KociolFwz(T,Twgx,Twzx,qK);
    Twgx=Twgn; Twzx=Twzn;
end
cicho=1;
% ...... Liczymy opoznienia 
Lrg=25; Lrz=15; dr=1.e-2; % dlug. i sredn rur w.g i w.z
sr=pi*dr^2/4; Vrg=Lrg*sr; Vrz=Lrz*sr; 
Fwmin=min(yFw);
dlg=Fwn/sr*dt/rog; dlz=Fwn/sr*dt/roz; dlmin=Fwmin/sr*dt/roz; 
drg=Vrg/Fwn*rog; ndrg=round(drg/dt);  if(ndrg<1) ndrg=1; end
drz=Vrz/Fwn*roz; ndrz=round(drz/dt);  if(ndrz<1) ndrz=1; end
Ndel=ceil(Lrg/dlmin); pozFwg(1:Ndel)=[1:Ndel]*dlg; pozFwz(1:Ndel)=[1:Ndel]*dlg; 
%return
% ...................................................
Tsym=12; Regul=0; 
ntSk=round(8/dt); DU=5; 
SymulGrzanie;
tx=t(ntSk-50: ntSk+1500);
Y=-(Tn(ntSk-50: ntSk+1500)-Tn(ntSk-50)); 
U=-(Twg(ntSk-50: ntSk+1500)-Twg(ntSk-50)); 
Ypocz=Tn(ntSk-1); Upocz=Twg(ntSk-1); Y=Y+Ypocz; U=U+Upocz; Unom=Upocz; 
nSkok=find(diff(U)>DU/2);
figure(111); plot(tx,Y,'k',tx,U*Kob+14.5,'b'); axis('tight'); xlabel('Tn(t) i Twg(t)');
Ldanych=length(Y); 
tSkok=nSkok*dt; stdV=0; rzad=3; DYshift=0; ntr=1; Unom=U(nSkok-1); Ypocz=Y(nSkok-1); 
Tokna=7*Taupom; Ysk0=Ypocz;
figure(111); plot(tx,Y,'k',tx,U+14.5,'b'); axis('tight'); xlabel('Tn(t) i Twg(t)');
tx=[0:Ldanych-1]*dt;
[hPom,Tos,Kos,dels,ds,Tor,Kor,delR,dr,Kob,DtR,nDt,nKpocz]=IdentOb(Ldanych,tx,dt,U,Y,Ysk0,tSkok,nSkok,Unom,Upocz,Ypocz,[],stdV,DU,DYshift,ntr,rzad,Kob,Taupom,Tokna);
delr=delR; 
% ..........................n.........................
Ld=length(Tzewn)-1; 
X(1)=T(1); X(2)=Ts(1); X(3)=Twz(1); X(4:lX)=Twg(1); 
X0=X; %Twg(1)=Twgn; qKa(1)=qKaln; Twzn(1)=Twz;
txthet='\Theta'; nrSym=1; txDel='\Delta';
dUtech=2.5; DUreg=0.05; Umin=25; Umax=90; Nielin=0;
Yref=ones(1,Ld+1)*T(1); DYref=-5; 
tS=round(Ld/4); TsRef=dt*tS; tS1=round(Ld*2/3); Yref(tS:tS1-1)=T(1)+DYref; Yref(tS1:Ld+1)=T(1); 
t=[0:Ld]*dt;
txPar=sprintf('rzad=%d K_{os}=%.2f/K_{ob}=%.2f T_{os}=%.2f/T_{ob}=%.2f d_s=%.2f (T_{os}+d_s)=%.2f',...
               rzad,Kos,Kob,Tos,Ts,dels,Tos+dels);
RegulPID;
vh=zeros(1,Ld+1); Nu=10; h=hPom; 
RegulPred;
