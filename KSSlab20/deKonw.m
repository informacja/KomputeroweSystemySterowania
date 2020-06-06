%deKonw.m
% ===== Identyfikacja nieparametryczna hp(i=1:Ldh))
clear hid Du y fif G FI Fi fiH yH Yob;
% ==============================================
stala=0; % to ma byc  0
bezh0=1; %0; %1; =0 z dobr¹ wag¹ dla h1 daje lepszy wynik
bezSel=1;
bezWag=0; wykl=0.6; %1.5; % jesli bezWag=0, to waga=(1/(Ldh-1))^(wykl/2)
Centr=1;
Hamm=1; %1; %1; %2;
if(Sigv<1.e-3) bezWag=1; %Hamm=0; 
else if(Sigv<2.01e-2) wykl=1;
    else if(Sigv<5.01e-2) wykl=0.6;
        else wykl=0.4;
        end
    end
end
tStud=1; if(tStud==1) lw11=4; else lw11=2; end
% ===============================================
Y0fig=Y0lin;
Ldr=length(Ur);
hpom=h; Ldh=length(hpom); %if(bezh0==0) Ldh=Ldh+1; end
if(jestDMC) Du=DUr; else Du=diff(Ur(1:Ldr)); end
Du=DUr; %Du=diff(Ur(1:Ldr));
Lid=length(Du);
if(Ldh>80) Ldh=80; end;
if(Ldr-2<3*Ldh)  Ldh=floor((Ldr-2)/3); end%if(Ldr-Ldh-2<2*Ldh)
%for(Ldhid=Ldh:-2:10)
%    if(Ldhid>Ldh) Ldhid=Ldh; end
    Nid=Ldr-Ldh-2;
    fif=zeros(Nid,Ldh);
    nid=[Ldh+1+bezh0:Ldr];
    if(bezSel) Nrid=[1:Nid];
    else
        % Yref(1:tS)=Y(1);Yref(tS+1:tS1)=1.4*Yref(tS); Yref(tS1+1:Ld+ntr)=0.6*Yref(tS);
        %Nrid=[Ldh+2+bezh0:100 200:Ldr-100]-Ldh-1-bezh0;
        Nrid=[Ldh+2+bezh0:80 200:Ldr-120]-Ldh-1-bezh0;
    end
    for(n=1:Nid)
        if(bezh0) %np=Ldh+n+1-bezh
            np=Ldh+n; tid(n)=(np+2)*DtR; % np wartoœæ n-1 dla y(n:Lid)
            %y0h(n,1)=Yr(np+1-Ldh); y(n,1)=Yr(np+2)-y0h(n,1);
            y(n,1)=Yr(np+1); %2);
            fif(n,1:Ldh-1)=Du(np:-1:np-Ldh+2); %:np-Ldh+1);
            %fif(n,Ldh)=Ur(np-Ldh+1); % +fif(n,end);
            fif(n,Ldh)=Ur(n+1); % +fif(n,end);
        else
            np=Ldh+n+1; tid(n)=(np+1)*DtR; % np wartoœæ n-1 dla y(n:Lid)
            %y0h(n,1)=Yr(np-Ldh); y(n,1)=Yr(np+1)-y0h(n,1);
            y(n,1)=Yr(np+1);
            fif(n,1:Ldh)=Du(np+1:-1:np-Ldh+2); %:np-Ldh);
            %fif(n,Ldh+1)=Ur(np-Ldh+1); %+fif(n,end);
            fif(n,Ldh+1)=Ur(n+2); %+fif(n,end);
        end
    end
    Ldh=Ldh+1-bezh0;
    figure(11+nf);
    subplot(lw11,1,1);
    th=[1:Ldh]*DtR; kol='brmgckbrmgckbrmgck';
    plot([2:nid(1)]*DtR,Yr([2:nid(1)])-Y0fig,'g',nid*DtR,Yr(nid)-Y0fig,'k',...
        [1:nid(1)-1]*DtR,Ur([1:nid(1)-1]),'c',(nid-1)*DtR,Ur(nid-1),'b'); %axis('tight');
    S=svd(fif'*fif); RCn=S(end)/S(1),
    % ......................................................
    clear FI Fi Yo ve vo;
    FI=fif; Yob=y; % FI i Yob - oryginalne wejœcia
    if(Centr)
        for(k=1:Ldh) fif(1:Nid,k)=fif(:,k)-mean(fif(:,k)); end;
        Yo=y-mean(y);
    else Yo=y;
    end
    % fif i Yo - scentrowane wejœcia i wyjscia
 dLdh=3; Ldhid=Ldh+dLdh;   
 %for(Ldhid=Ldh:-3:10)  
 powt=0;
 while(Ldhid>10)
    clear hid G Fi fiH yH JHmin JH; 
    Ldhid=Ldhid-dLdh; 
    JHM=1.e20;
    if(Ldhid>Ldh) Ldhid=Ldh; end
    Fi=fif(:,1:Ldhid); % ograniczamy identyfikacjê do h(1:Ldhid) punktow
    % Yo - scentrowane wyjœcia; Fi scentrowane wejœcia dla h(1:Ldhid)
    if(Ldhid<Ldh)
        % przeliczenie Ldhid na bez h(0):
        %                     bezh=0 (jest h(0)) Ldhf=Ldhid-1
        %                     bezh=1 (zaczynamy od h(1)) Ldhf=Ldhid
        for(n=1:Nid)
            %if(bezh0) % zaczynamy od h(1)
            %    np=Ldh+n; tid(n)=(np+2)*DtR;
            %    fif(n,Ldh)=Ur(np-Ldhf+1); % +fif(n,end);
            %else
            %np=Ldh+1-bezh+n
            np=Ldh+n; tid(n)=(np+1)*DtR; % np wartoœæ n-1 dla y(n:Lid)
            %Fi(n,Ldhid)=Ur(n+1);
            Fi(n,Ldhid)=Ur(np-Ldhid);
            %end
        end
        if(Centr) Fi(:,Ldhid)=Fi(:,Ldhid)-mean(Fi(:,Ldhid)); end
    end
    % ........................................................
    clear Rh Rh1;
    if(bezWag) Rh=zeros(Ldhid+stala,Ldhid+stala);
    else
        R1h=zeros(Ldhid,Ldhid);
        for(i=1:Ldhid-1)
            %r=(1/(Ldh-i)); %
            r=(1/(Ldhid-i))^(wykl/2); % najlepszy wariant
            %r=sqrt(i/(Ldh-1));
            R1h(i,i)=r; R1h(i,i+1)=-r;
        end
        %if(bezh0==0) R1h(1,1)=R1h(1,1)+1; else R1h(1,1)=1; R1h(1,2)=0; end
        R1h(1,1)=R1h(1,1)+0.1;
        Rh=R1h'*R1h;
        if(stala) Rh(:,Ldhid+1)=0; Rh(Ldhid+1,:)=0; end
    end
    % ....................................................
    Niden=Nid;
    if(Hamm==0)
        fi=Fi(Nrid,:); y=Yo(Nrid,1); Nid=length(y(:,1));
        fiW=fi; litR=5; if(Sigv<1.e-3) litR=1; end
    else % model Hammersteina
        litR=0;
        litH=50; a10=0.7; a1p=a10+0.04; a1f=0.99; da1=(a1f-a10)/(litH-1); a1=a10;
        a20=0; a2p=0.04; a2f=0.99; da2=(a2f-a2p)/(litH-1); a2=a2p;
        if(Hamm>1) litH2=litH; else litH2=1; end
        figure(234+Centr+bezSel); Kol='cgbmrk'; clear JH;
        a1=a10;   JHM=1.e20;
        for(it=1:litH+1)
            clear yH fiH fi y; 
            a2=a20; %a2p;
            for(m=1:litH2)
                yH=Yo(1+Hamm:Nid,1)-Yo(Hamm:Nid-1,1)*a1;
                fiH(1:Nid-Hamm,:)=Fi(1+Hamm:Nid,1:Ldhid)-Fi(Hamm:Nid-1,1:Ldhid)*a1;
                lfiH=length(fiH(:,1));
                if(Hamm>1)
                    yH=yH-Yo(1:Nid-2,1)*a2;
                    fiH(1:Nid-2,:)=fiH(1:Nid-2,:)-Fi(1:Nid-2,1:Ldhid)*a2;
                end
                fi=fiH(Nrid(1:lfiH),:); y=yH(Nrid(1:lfiH),1);
                if(stala)
                    fi(1:length(fi(:,1)),Ldhid+1)=1;
                end
                G=fi'*fi;
                if(1) hid=inv(G+Rh)*fi'*y; else hid=inv(G)*fi'*y; end
                Yh=fi*hid;
                E=y-Yh; lE=length(E);
                JH(it,m)=(E'*E)/lE;
                if(JH(it,m)<JHM)
                    if(stala) v0=hid(Ldhid+1); hm=hid(1:Ldhid); clear hid; hid=hm; end
                    hm=hid; A1=a1; A2=a2; JHM=JH(it,m); JHmin(it)=JHM; srE=mean(E),;
                    if(it<7) kol=[Kol(it) '.-'];
                    else if(it<13) kol=[Kol(it-6) '.:'];
                        else if(it<19) kol=[Kol(it-12) '.--'];
                            else kol=[Kol(mod(it-1,6)+1) '.-.'];
                            end
                        end
                    end
                    subplot(3,1,1); Ldhp=min(Ldhid,length(hpom));
                    plot([1:Ldhid]*DtR,hid,kol,[1:Ldhp]*DtR,hpom(1:Ldhp),'r*-');
                    xlabel(sprintf('a1=%.2f a2=%.2f',A1,A2)); hold on;
                end
                a2=a20+(m-1)*da2;
            end
            %if(it>6) i=mod(it-1,6)+1; kol=[Kol(i) '.:']; else kol=[Kol(it) '.-']; end
            if(it<7) kol=[Kol(it) '.-'];
            else if(it<13) kol=[Kol(it-6) '.:'];
                else if(it<19) kol=[Kol(it-12) '.--'];
                    else kol=[Kol(mod(it-1,6)+1) '.-.'];
                    end
                end
            end
            if(Hamm>1)
                subplot(3,1,2); plot(JH(it,:),kol); hold on;
            end
            %nH=find(JH(it,:)==min(JH(it,:))); A1=a1; A2=(nH-1)*da2; hid=hm;
            %JHmin(it)=JH(it,nH); A1min(it)=A1; A2min(it)=A2;
            %subplot(3,1,1); plot([1:Ldh]*DtR,hid,'k.-',[1:Ldh]*DtR,hpom(1:Ldh),'r*-');
            LYo=length(Yo); tid=nid(3-bezh0:end)*DtR;
            ve=Yo-Fi*hm;  vo(1,1)=ve(1,1); vo(2:LYo,1)=ve(1:LYo-1,1)*A1;
            subplot(3,1,3); hold off; plot(tid,ve,'k',tid,vo,'r:',tid,ve-vo,'g');
            xlabel(sprintf('a1=%.2f a2=%.2f A1=%.2f A2=%.2f',a1,a2,A1,A2)); axis('tight'); hold on;
            a1=a1+da1; %0+it*da1; %+da1; a2=a2p;
        end
        hid=hm; hid(Ldhid+1:Ldh)=hm(Ldhid);
        a1=A1; a2=A2; 
        %nHm=find(JHmin==min(JHmin)); A1m=A1min(nHm); A2m=A2min(nHm);
        Ldhp=min(Ldh,length(hpom));
        subplot(3,1,1); hold off; plot([1:Ldh]*DtR,hid,'k.-',[1:Ldhp]*DtR,hpom(1:Ldhp),'r*-');
        xlabel(sprintf('a1=%.2f a2=%.2f',A1,A2)); hold off;
        subplot(3,1,2);
        if(Hamm==1) hold off; plot(a10+[0 1:litH]*da1, JH,'k.-'); end
        axis('tight'); xlabel(sprintf('J(a1,a2): JHmin=%.3g a1m=%.2f a2m=%.2f',JHM, A1, A2)); hold off;
        subplot(3,1,1); axis('tight'); xlabel(sprintf('h(a1,a2): a1m=%.2f a2m=%.2f',A1, A2)); hold off;
        Nid=length(y(:,1));
    end
    figure(11+nf);
    subplot(lw11,1,1); hold off;
    clear hm;
    for(it=1:litR) % jesli model Hammersteina litR=0
        %            obliczenia wg uogóln.MNK
        G=fiW'*fi;
        if(1) hm=inv(G+Rh)*fiW'*y; else hm=inv(G)*fiW'*y; end
        Yh=fi*hm;
        Eid=y-Yh; le=length(Eid);
        hid=hm; hid(Ldhid+1:Ldh)=hm(Ldhid);
        Yhid=FI*hid;
        E=Yob-Yhid; Esr=mean(E); Ec=E-Esr;
        lE=length(E); alfE=(Ec(2:lE)'*Ec(1:lE-1))/(lE*var(Ec));
        RE=eye(lE); for(i=1:lE) for(j=i+1:le) RE(i,j)=alfE^(j-i); RE(j,i)=RE(i,j); end, end
        W=inv(RE); fiW=W'*fi;
        %e=E(2:lE)-RE*E(1:lE-1);
        %figure(301); plot([1:lE],E,'k',[2:lE],e,'r'); axis('tight')
        %le=length(e); esr=mean(e),; ec=e-esr;
        %alfe=(ec(2:le)'*ec(1:le-1))/(le*var(ec))
        subplot(lw11,1,1);
        plot(nid(3-bezh0:end)*DtR,Yhid,[kol(it) ':'],nid(3-bezh0:end)*DtR,E,'b:');
        hold on;
        subplot(2,1,2); %hold on;
        if(bezh0)
            plot([DtR th],[0; hid],kol(it)),
        else plot(th,hid,kol(it)),%th,[0 hpom(1:Ldh-1)],'r.-');
        end
        hold on;
    end
    Koniec=1; if(Ldh-Ldhid<2*dLdh) Koniec=0; continue; end
    for(i=Ldhid-1:-1:1) if((hid(i)/hid(Ldhid))>1.01) Koniec=0; break; end, end
    if(Koniec) 
        npl=find(hid(1:Ldhid)/hid(Ldhid)>0.5); if(isempty(npl)) np=2; else np=npl(1); end
        if(np<2) np=2; end
        dhw=diff(hid(np:Ldhid))/hid(Ldhid); npl=find(dhw<-0.1); 
        if(isempty(npl)) 
            if(dhw(end)>dhw(end-1) && powt) Ldhid=Ldhid+round(dLdh/2); Koniec=0; powt=1; else break; end
        else Koniec=0; end
    end
end
subplot(2,1,2); hold off;
th=[1:Ldh]*DtR; 
if(tStud) 
% ............. Obliczamy statyst. Studenta ..................
    clear tSt Ka sigh;
    lF=length(Fi(1,:)); 
    Ka=inv(G)*(E'*E)/(lE-lF); 
    for(i=1:lF) sigh(i,1)=sqrt(Ka(i,i)); tSt(i)=abs(hid(i))/sigh(i,1); end
    sigh(Ldhid+1:Ldh,1)=sigh(Ldhid,1); tSt(Ldhid+1:Ldh)=tSt(Ldhid); 
    subplot(lw11,1,2); 
    if(lw11==2)
        hold on; plot(th,hid,'k.-',th,hid+sigh,'c:',th,hid-sigh,'c:'); axis('tight'); hold off; 
    else
        nh=find(tSt>1); if(isempty(nh)) nis=1; else nis=nh(1); end        
        subplot(lw11,1,2); plot(th,tSt,'k.-'); axis('tight'); 
        ax=axis; hold on; plot([nis nis]*DtR,ax(3:4),'r--'); hold off; 
        xlabel('Statystyki Studenta dla odpowiedzi skokowej');
    end
end
subplot(lw11,1,1); hold on;
Yhid=FI*hid;
if(bezh0) hid=[0 hid']; th=[0 1:Ldh]*DtR; Ldh=Ldh+1; else hid=hid'; end
plot([2:nid(1)]*DtR,Yr([2:nid(1)])-Y0fig,'g',[1:nid(1)-1]*DtR,Ur([1:nid(1)-1]),'c',...
    nid(3-bezh0:end)*DtR,Yhid,'r',nid*DtR,Yr(nid)-Y0fig,'k',...
    (nid-1)*DtR,Ur(nid-1),'b'); %axis('tight');
axis('tight');
%ax=axis; hold on; plot(nid(2:end)*DtR,y0h,'m:',nid([1 1])*DtR,ax(3:4),'r--'); hold off;
nEmax=find(max(abs(E))==abs(E))+1;
ax=axis; hold on; plot(nid([1 1])*DtR,ax(3:4),'r--',nid([nEmax nEmax])*DtR,ax(3:4),'g:'); hold off;
xlabel(sprintf('Dane z %s dla odpow.skokowej: L_{dr}=%d N_{id}=%d L_{dh}=%d',plikY,Ldr,Nid,Ldh));
subplot(2,1,2); if(litR>1) hold on; else hold off; end
Ldhp=min(Ldh,length(hpom));
%plot(th,hid,'k.-',[0 1:Ldhp]*DtR,[0 hpom(1:Ldhp)],'r.-');
plot(th,hid,'k.-',[0 1:Ldhp-1]*DtR,[hpom(1:Ldhp)],'ro-');
if(tStud>1) hold on; plot(th,hid-sigh,'m.:',th,hid+sigh,'m.:'); end
axis('tight'); ax=axis; hold on; plot(ax(1:2),[0 0],'k:'); hold off;
xlabel(sprintf('Regresyjna (k.-) i zmierzona (r.-) odpowiedz skokowa: N_{id}=%d L_{dh}=%d',Nid,Ldh));
%subplot(3,1,3); plot(th,hpom(1:Ldh),'k.-'); axis('tight'); xlabel(sprintf('Zmierzona odpowiedz skokowa: L_{dh}=%d',Ldh));
npp=find(hid>0); if(isempty(npp)) np=1; else np=npp(1); end
nhk=min(length(hid),length(hpom)); np=1'
Eh=hid(np:nhk)-hpom(np:nhk); sigEh=std(Eh);
%idenh=0;
if(sigEh<sEhMin)
    sEhMin=sigEh; idenh=1; nwar(idenh)=nw;  Nfig(idenh)=nf;
    mwar(idenh)=m;  Yhiden(idenh,:)=hid;
else
    if(sigEh<sehMin)
        sehMin=sigEh; idenh=2; Nfig(idenh)=nf;
        nwar(idenh)=nw; mwar(idenh)=m; Yhiden(idenh,:)=hid;
    end
end
return
