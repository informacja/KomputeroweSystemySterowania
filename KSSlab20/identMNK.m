% plik identMNK.m
%plikU='Uo_regPID.csv'; plikY='Yo_regPID.csv'; plikYref='Yref_regPID.csv';
%U=csvread(plikU); Y=csvread(plikY); Yref=csvread(plikY);
%xtu=[1:length(U)]*dt; ntY=[1:length(Y)]; xty=ntY*dt; plot(xtu,U,'b',xty,Yid,'k',xty,Yref(ntY),'r'); axis('tight');
% ------------- Szukanie optymalnego okna -----------------
global przyr wagi nfig1; 
if(skok)
    d1=2; 
    Ld=Ldanych;
    Ldr=round(Ld*0.25);  %Ldr=round(Ld*0.2);
    fig100=0;
    if(fig100) figure(100); end
    clear RCN Kyu Kyy;
    RCNmax=-1; dmax=d1; koniec=0;
    for(d=d1:Ld)
        ldf=Ldr+d-1; if(ldf>Ld) break; end
        clear y fi A G okno;
        okno=[d:ldf];
        %y=Y(okno)'; fi(:,1)=Y(okno-1)'; fi(:,2)=U(okno-1)';
        y0=Y(okno(1)-1); u0=U(okno(1)-1);
        y=Y(okno)'-y0; fi(:,1)=Y(okno-1)'-y0; fi(:,2)=U(okno-1)'-u0;
        G=fi'*fi;
        S=svd(G); RCN(d-1)=S(end)/S(1); dRCN=RCN(d-1); %-RCN(d1-1);
        ysr=mean(y); usr=mean(fi(:,2)); ypsr=mean(fi(:,1));
        
        Kyu(d-1)=((y-ysr)'*(fi(:,2)-usr))/(length(y));
        Kyy(d-1)=((y-ysr)'*(fi(:,1)-ypsr))/(length(y));
        if(dRCN>RCNmax) RCNmax=dRCN; dmax=d; end
        if(d==nSkok-1 && skok) dpo=d; koniec=1; end
        %    else
        %        if(dRCN<RCNmax*0.5 && d-dmax>100)
        %            if(okno(end)*dt>1.5*tSkok) dpo=d; koniec=1; end
        %        end
        %    end
        if(fig100&&(mod(d-2,100)==0 || koniec))
            subplot(2,1,1);
            plot(okno*dt,y,'k',(okno-1)*dt,fi(:,2),'r'); axis('tight'); hold on;
            xlabel(sprintf('Okno obserw. od t=%.2f do %.2f; y0=%.2f u0=%.2f',okno(1)*dt, okno(end)*dt,y0,u0));
            ax=axis; plot(okno([1 1])*dt,ax(3:4),'r:',dt*okno([Ldr Ldr]),ax(3:4),'k:'); axis([0.999*ax(1) 1.001*ax(2) ax(3:4)]);
            hold off;
            xd=[d1:d-1]*dt; subplot(2,1,2); plot(xd,RCN([d1:d-1]),'k');
            if(koniec) dpo=d; break; end
            %input(' ?? ');
        end
        if(koniec) dpo=d; break; end
    end
    if(fig100)
        hold off;
        xd=[d1:d-1]*dt;
        lw=4; lk=1; lfp=1;
        subplot(lw,lk,1);
        plot(okno*dt,y,'k',(okno-1)*dt,fi(:,2),'r'); axis('tight'); hold on;
        xlabel(sprintf('Okno obserw. od t=%.2f do %.2f; y0=%.2f u0=%.2f',okno(1)*dt, okno(end)*dt,y0,u0));
        ax=axis; plot(okno([1 1])*dt,ax(3:4),'r:',dt*okno([Ldr Ldr]),ax(3:4),'k:'); axis([0.999*ax(1) 1.001*ax(2) ax(3:4)]);
        hold off;
        subplot(lw,lk,lfp+1); plot(xd,RCN(d1-1:d-2),'k');
        axis('tight'); xlabel(sprintf('RCN dla okno=%.2f',Ldr*dt));
        subplot(lw,lk,lfp+2); plot(xd,Kyu(d1-1:d-2),'k');
        axis('tight'); xlabel(sprintf('Kyu dla okno=%.2f',Ldr*dt));
        subplot(lw,lk,lfp+3); plot(xd,Kyy(d1-1:d-2),'k');
        axis('tight'); xlabel(sprintf('Kyy dla okno=%.2f',Ldr*dt));
    end
    if(skok==0); nd=find(RCN/min(RCN)>7); dpo=nd(1); if(dpo<d1) dpo=d1;end, end
    %return;
    d1=dpo; td1R=d1*dt;   
    Y0=mean(Y(1:dpo)); U0=mean(U(1:dpo));
    dMax=(Ld-Ldr)/3+d1; 
    if(dMax>Ldr/2+d1) dMax=floor(Ldr/2)+d1; end
    if(dMax+Ldr>Ld) dMax=Ld-Ldr; end
    Ldmod=Ldanych; nDtR=nDt; 
    %Ldy=floor(Ldr/nDt); W=eye(Ldy);
    %if(wagi) for(i=1:Ldy) W(i,i)=(1/(Ldy+1-i))^4; end; end
%else
end
nskok=2; % pozycja skoku w tablicy 
      Kor=Kob; Tor=Tos; delR=dels; dr=ds; delr=dels; 
      d01=1; dmin=d01; delmin=d01-2*nDt; if(delmin<1) delmin=1; end
      Ndrf=Ldyh+delmin; 
      [Kor,Tor,dr,delr,Ndrf]=idhRegr(yh1(nSkok-(nskok-1)*nDt:end),U(nSkok-(nskok-1)*nDt:end)-u0,0,0,dt,nDt,DtR,1,dMax,Ldyh,Ld,nskok,Kos);
      % Tu wywolujemy yh1 1.krok przed skokiem. Skok jest w chwili nSkok
      % a yh1(nSkok-ndT+nDt) jest wart. w chwili skoku    
% ...... Koniec identyfikacji regres. ........
%Ypr=Unom*Kor; Yps=Unom*Kos;
if(dr<=0) dr=1; end; if(ds<=0) ds=1; end
Ypr=0; Yps=0; % Obliczenia bêd¹ dla odcy³ki od wart.nominalnej
clear Yr Ys t;
% === Przebieg symulacji dla narysowania wynikow identyfikacji 
t(1:ntr)=[0:ntr-1]*dt; Yr(1)=Ypr; Ys(1)=Yps; 
for(n=1:Ldmod) 
    t(n+ntr)=(n+ntr)*dt;
    Yr(n+ntr)=model(Kor,Tor,dr,U-Unom,n,Ypr); Ypr=Yr(n+ntr);  
    Ys(n+ntr)=model(Kos,Tos,ds,U-Unom,n,Yps); Yps=Ys(n+ntr);   
end
Yr=Yr+Ypocz; Ys=Ys+Ypocz; % dodanie wartoœci nominalnej
%Ys=Ys-Ysk0; Yr=Yr-Ysk0; % przesuniêcie dla lepszej wizualizacji
figure(nfig1);
subplot(2,1,1); 
hold on; plot(t,Ys-Ysk0,'r'); %,t,Yr,'m'); 
dhmin=Ypocz-Ysk0; 
if(skok)
    %Ysk0=Ypocz-Y0lin;
    ax=axis; 
    plot([tSkok tSkok],ax(3:4),'b:',[npocz npocz]*dt,ax(3:4),'m:',...
        [nkonc nkonc]*dt,ax(3:4),'m:',[nKpocz nKpocz]*dt,ax(3:4),'b--'); 
    plot([npocz nkonc]*dt,[hmin hmin]+dhmin,'m:',[npocz nkonc]*dt,[hmax hmax]+dhmin,'m:')
    txPar=sprintf('rzad=%d K_{os}=%.2f/K_{ob}=%.2f T_{os}=%.2f/T_{ob}=%.2f d_s=%.2f (T_{os}+d_s)=%.2f',...
               rzad,Kos,Kob,Tos,Ts,dels,Tos+dels);
    xlabel(txPar); 
    hold off; 
else n=1; 
end
subplot(2,1,2); hold off; 
nU=1:length(U); 
if(~isempty(Yt))
    plot(t(nU),U(nU),'g',t(nU),Y(nU)-Ysk0,'k',t(nU),Yt(nU)-Ysk0,'b--',t,Yr(nU)-Ysk0,'r',[1:length(yh1)]*dt,yh1+Ypocz-Ysk0,'r:'); 
else plot(t(nU),(U(nU)-Upocz)*Kor,'g',t(nU),Y(nU)-Ysk0,'k',t(nU),Yr(nU)-Ysk0,'r',[1:length(yh1)]*dt,yh1+Ypocz-Ysk0,'r:'); 
end
axis('tight'); ax=axis; hold on; 
tkR=td1R+Ldr*dt+delr;
plot([td1R td1R],ax(3:4),'b:',[td1R td1R]+delr,ax(3:4),'r:',[tkR tkR],ax(3:4),'r:',tSkok+[Ndrf Ndrf]*dt,ax(3:4),'k:'); 
txPar=sprintf('rzad=%d K_{or}=%.2f/K_{ob}=%.2f T_{or}=%.2f/T_{ob}=%.2f d_r=%.2f (T+d)=%.2f',...
           rzad,Kor,Kob,Tor,Ts,delr,Tor+delr); 
xlabel(txPar); axis('tight');
hold off; 
%return
