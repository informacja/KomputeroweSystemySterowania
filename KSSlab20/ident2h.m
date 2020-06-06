%ident2h
if(idenh>0) idM=2; else idM=0; end
for(id=1:idM) % Identyfikacja odpow. skokowej
    idenh=id; clear Y he;
    Y=Yhiden(idenh,:);
    Uh=[0 1:Ldh-1];
    Ldmod=length(Uh);
    % Identyfikacja modelu met. aproksymacji
    m=0; npocz=0; nkonc=Ldmod; DU=1;
    clear y he;
    Ksk(idenh)=Y(end); hMax=-1.e20;
    %Ks=max(Y); Ks=Ks(1); he=Y*Ksk(idenh)/Ks;
    %np=find(he>0); nplus=np(1);
    np=find(Y>Ksk(idenh)); 
    if(isempty(np)) n1=length(Y); else n1=np(1); end
    Ks=mean(Y(n1:end)); 
    he=Y(1:n1-1); he(n1:Ldmod)=hsr; 
    he=he*Ksk(idenh)/Ks; 
    for(n=1:Ldmod)
        hY=he(n);
        if(hY>hMax && m==0) hMax=hY; nmax=n; end
        if(hY>0.1*Ksk(idenh)) %* && hY<hmax)
            if(hY<hMax) m=0; npocz=n+1;
            else
                if(hY>0.98*Ksk(idenh))
                    nkonc=n-1; break;
                end
                m=m+1;
                y(m)=log(1-hY/Ksk(idenh)); % obserwacje wyjúcia
            end
        else
            if(m>0)
                m=0; clear y;
            else
                if(m==0) npocz=n+1;
                else
                    nkonc=n-1;
                end,
            end
        end
    end
    figure(113); % nkonc=Ldmod; %nkonc+1;
    subplot(2,1,1);
    plot([1:npocz]*DtR,he(1:npocz),'go-',[npocz:nkonc]*DtR,he(npocz:nkonc),'ro-',...
        [1:npocz]*DtR,Y(1:npocz),'g:',[npocz:Ldmod]*DtR,Y(npocz:Ldmod),'r:');
    if(npocz==0) npocz=1; fprintf(1,'\nZa silne zaklocenie !!! npocz=%d',npocz); end
    % Wyznaczamy sta≥π czasowπ i opoünienie regresjπ liniowπ
    tm=[1:m]*DtR; ty=[npocz:nkonc]*DtR;
    %if(nkonc<Ldmod) y=[y log(1-he(Ldmod)/Ksk(idenh))]; tm=[tm Ldmod*DtR]; ty=[ty Ldmod*DtR]; end
    subplot(2,1,2);plot(ty,y,'ro-');
    a=polyfit(ty,y,1); % a(1)=-1/To; a(2)=d/To;
    subplot(2,1,2); plot(ty,y,'ko-',ty,a(1)*ty+a(2),'r');
    Tsk(idenh)=-1/a(1); %*Ks/Ksk(idenh); 
    dels=a(2)*Tsk(idenh); %-npocz*DtR;
    if(dels<0) dels=DtR; end
    dsk(idenh)=round(dels/DtR); 
    ds=dsk(idenh); Kos=Ksk(idenh); Tos=Tsk(idenh);
    delS=ds*DtR; dSk(idenh)=delS; 
    clear Uh;
    Uh(1)=0; Uh(2:Ldmod)=1; d1=1; dMax=ds+3*round(Tos/DtR); 
    Ldid=Ldmod-dMax-1; 
    Ldr=Ldhid; %+round((Ldmod-Ldhid)/5); 
    if(Ldr<10) Ldr=10; dMax=Ldmod-Ldr-1; end
    jestK=Y(end);
    %[Kor,Tor,dr,delr,nd]=idhRegr(Y,Uh,y0,u0,dt,nDt,DtR,d1,dMax,Ldr,Ld,nSkok,jestK)
    [Kor,Tor,dr,delR,Ndrf]=idhRegr([0 Y],Uh,0,0,DtR,1,DtR,1,dMax,Ldr,Ldmod,2,jestK); 
    dreg(idenh)=dr; Kreg(idenh)=Kor; Treg(idenh)=Tor;
    %delR=dr*DtR; 
    dReg(idenh)=delR; 
    clear Yrg Ysk t;
    % === Przebieg symulacji dla narysowania wynikow identyfikacji
    %t(1:ds)=[0:ds-1]*DtR; 
    ndod=1;
    Ysk(1:ndod)=0; Yps=0; Ypr=0; Yrg(1:ndod)=0; 
    Ut=[0 ones(1,Ldmod*nDt-1)];
    ndod=1; t(1)=0; dtim=dt; %dt=DtR; 
    for(n=1:Ldmod*nDt)
        nt=n+ndod;
        t(nt)=(nt)*DtR;
        Yrg(nt)=model(Kor,Tor,dr*nDt,Ut,n,Ypr); Ypr=Yrg(nt);
        Ysk(nt)=model(Kos,Tos,ds*nDt,Ut,n,Yps); Yps=Ysk(nt);
    end
    dt=dtim;
    subplot(2,1,1); hold on;
    plot([1:nDt:length(Ysk)]*dt,Ysk(1:nDt:end),'b.--'); 
    plot([1:nDt:length(Yrg)]*dt,Yrg(1:nDt:end),'r.--'); 
    hold off; axis('tight'); hold off;
    figure(11+Nfig(idenh));
    subplot(2,1,2); hold on;
    %plot([0:length(Ys)-1]*DtR,Ys,'m+-',[0:length(Yr)-1]*DtR,Yr,'mx-',[0:Ldmod-1]*DtR,Y,'k.-'); 
    plot([2:length(Ysk)]*dt,Ysk(2:end),'m',[2:length(Yrg)]*dt,Yrg(2:end),'m',[0:Ldmod-1]*DtR,Y,'k.-'); 
    hold on;
    plot([1:nDt:length(Ysk)]*dt,Ysk(1:nDt:end),'m+',[1:nDt:length(Yrg)]*dt,Yrg(1:nDt:end),'mx'); 
    axis('tight'); ax=axis; plot([Ndrf Ndrf]*DtR,ax(3:4),'k:');
    hold off; txth='\theta';
    xlabel(sprintf('Pom ro-; Obl k.-; Apr m+-; Reg mx; odp.skok.: N_{id}=%d L_{dh}=%d; K_{os}=%.3f T_{os}=%.3f %s_{os}=%.3f; K_{or}=%.3f T_{or}=%.3f %s_{or}=%.3f',Nid,Ldh,Kos,Tos,txth,delS,Kor,Tor,txth,delR));
end
