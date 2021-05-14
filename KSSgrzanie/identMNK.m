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
    figure(100);  clear RCN Kyu Kyy;
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
        if(mod(d-2,100)==0 || koniec)
            subplot(2,1,1);
            plot(okno*dt,y,'k',(okno-1)*dt,fi(:,2),'r'); axis('tight'); hold on;
            xlabel(sprintf('Okno obserw. od t=%.2f do %.2f; y0=%.2f u0=%.2f',okno(1)*dt, okno(end)*dt,y0,u0));
            ax=axis; plot(okno([1 1])*dt,ax(3:4),'r:',dt*okno([Ldr Ldr]),ax(3:4),'k:'); axis([0.999*ax(1) 1.001*ax(2) ax(3:4)]);
            hold off;
            xd=[d1:d-1]*dt; subplot(2,1,2); plot(xd,RCN([d1:d-1]),'k');
            if(koniec) dpo=d; break; end
            %input(' ?? ');
        end
    end
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
%if(1) [Kor,Tor,dr,delr,nd]=idhRegr(Y,U(,y0,u0,dt,nDt,DtR,d1,dMax,Ldr,Ld,nSkok,Kos);
nskok=2; % pozycja skoku w tablicy 
if(1) Kor=Kob; Tor=Tos; delR=dels; dr=ds; delr=dels; 
      d01=1; dmin=d01; delmin=d01-2*nDt; if(delmin<1) delmin=1; end
      Ndrf=Ldyh+delmin; 
      [Kor,Tor,dr,delr,Ndrf]=idhRegr(yh1(nSkok-(nskok-1)*nDt:end),U(nSkok-(nskok-1)*nDt:end)-u0,0,0,dt,nDt,DtR,1,dMax,Ldyh,Ld,nskok,Kos);
      % Tu wywolujemy yh1 1.krok przed skokiem. Skok jest w chwili nSkok
      % a yh1(nSkok-ndT+nDt) jest wart. w chwili skoku    
else
    clear RCN J Kyu Kyy;
    nd=0; R=[2 0; 0 1]; R=[1 0; 0 1]*0.;
    Jdmin=1.e20; Jdmax=0; RCNmax=-1;
    for(d=d1:dMax)
        nd=nd+1; %nd=d-d1+1;
        dyf=d+Ldr; if(dyf>Ld) dyf=Ld; break; end
        nuf=d1+Ldr-1; nuf=d1+Ldr;
        clear y fi A G yn yp ud;
        %y0=Y(d); u0=U(d1);
        yn=Y(d+nDt:nDt:dyf)'-y0;
        yp=Y(d:nDt:dyf-nDt)'-y0; ud=U(d1:nDt:nuf-nDt)'-u0;
        if(przyr)
            y=diff(yn); yn=y; % Y(d+1:d+Ldr+1-1)
            fi(:,1)=diff(yp); fi(:,2)=diff(ud);
        else
            y=yn; fi(:,1)=yp; fi(:,2)=ud;
        end
        clear yp ud;
        Fi=fi;
        kt=1; deltR=floor(nd/nDt);
        km=round(deltR/4); if(km==0) km=1; end
        for(k=1:km)
            y=yn(k:end,1); fi=Fi(k:end,:);
            lfi=length(fi(:,2)); W=eye(lfi);
            if(0) %k>1)
                fi(:,1)=yt(1:lfi);
            end
            if(wagi==0) fiW=fi'; GW=fi'*fi;
            else
                if(k>1)
                    alf=A(1); Kor=A(2)/(1-alf); Tor=-DtR/log(alf);
                    nWag=round(Tor/DtR);
                else nWag=deltR;
                end
                nWag=min(nWag,lfi);
                for(i=1:nWag) W(i,i)=(1/(nWag+1-i))^4; end;
                fiW=fi'*W; GW=fiW*fi;
            end
            S=svd(GW); RCNx(k)=S(end)/S(1);
            A=inv(GW+R)*fiW*y; %if(A(1)>1) continue; end
            yp=fi(1,1); clear yt; yt=fi*A;
            %for(i=1:lfi) yt(i,1)=yp*A(1)+fi(i,2)*A(2); yp=yt(i,1); end
            e=y-yt;
            Jx(k)=(e'*e)/length(yt);
            ysr=mean(y); usr=mean(fi(:,2)); ypsr=mean(fi(:,1));
            %Kyux(k)=((y-ysr)'*(fi(:,2)-usr))/(length(y));
        end
        shift=k-1;
        delr=nd*dt; alf=A(1); Kor=A(2)/(1-alf); Tor=-DtR/log(alf);
        %clear yn Fi;
        J(nd)=(e'*W*e)/length(yt);
        RCN(nd)=S(end)/S(1);
        if(RCN(nd)>RCNmax) RCNmax=RCN(nd); dRCmax=nd; end
        if(przyr)
            sgy=std(y); sgu=std(fi(:,2)); sgyp=std(fi(:,1));
        else sgy=1; sgu=1; sgyp=1;
        end
        %Kyu(nd)=((y-ysr)'*(fi(:,2)-usr))/(length(y)*sgy*sgu);
        Kyy(nd)=((y-ysr)'*(fi(:,1)-ypsr))/(length(y)*sgy*sgyp);
        if(J(nd)>Jdmax) Jdmax=J(nd); dmax=nd; Jdmin=1.e20;
        else
            if(J(nd)<Jdmin)
                dmin=nd; Jdmin=J(nd);
                Ao=A; fio=fi; Fio=Fi; yno=yn; yo=y;
            else
                if(Jdmin>0 && delr>1.5*Tor) break; end
                if(0) %nd-dmin>10 && dmax<dmin)
                    if(Jdmax/Jdmin>1.50) break; end
                end
            end
        end
        if(mod(nd,100)==0)
            nsh=shift*nDt;
            tSk=ceil((nSkok-d1+1)/nDt)*DtR; %tSko=tSk-d+1;
            figure(200); subplot(3,1,1);
            xdy=[1: length(yn)]*DtR; xdu=xdy;
            
            LFi=length(Fi(:,1)); yp=Fi(1,1);
            for(i=1:LFi) ytf(i)=yp*A(1)+Fi(i,2)*A(2); yp=ytf(i); end
            plot(xdy,yn,'k.-',xdu,Fi(:,2),'r.-',xdy,ytf,'m.'); axis('tight');
            ax=axis; hold on;
            plot([tSk tSk], ax(3:4),'b:',[tSk tSk]+nsh, ax(3:4),'b-.'); hold off;
            axis([0.99*ax(1) 1.01*ax(2) ax(3:4)]);
            
            subplot(3,1,2);
            %xdy=[d+(przyr+1+shift)*nDt:nDt:dyf]*dt; xdu=[d1+(przyr+shift)*nDt:nDt:nuf-nDt]*dt;
            xdy=[d+(przyr+1)*nDt:nDt:dyf]*dt; xdu=[d1+(przyr)*nDt:nDt:nuf-nDt]*dt;
            ny0=[d1+(przyr)*nDt:nDt:d+(przyr+1)*nDt];
            
            yp=fi(1,1); clear ytm;
            for(i=1:lfi) ytm(i,1)=yp*A(1)+fi(i,2)*A(2); yp=ytm(i,1); end
            
            plot(ny0*dt,Y(ny0)-y0,'g.-',xdy(shift+1:end),yt,'y.:',xdy(shift+1:end),ytm,'c.-',xdy,yn,'k.-',xdu,Fi(:,2),'r.-',xdy,ytf,'m.'); axis('tight');
            hold on; ax=axis;
            plot([tSkok tSkok], ax(3:4),'b--',[d+nDt d+nDt]*dt,ax(3:4),'g:',[d+nDt+nsh d+nDt+nsh]*dt,ax(3:4),'b-.',...
                [dyf dyf]*dt,ax(3:4),'b:',[d1 d1]*dt,ax(3:4),'r:',[nuf-nDt nuf-nDt]*dt,ax(3:4),'r:');
            axis([0.99*ax(1) 1.01*ax(2) ax(3:4)]); hold off;
            delr=nd*dt; alf=A(1); Kor=A(2)/(1-alf); Tor=-DtR/log(alf);
            xlabel(sprintf('nd=%d J=%.3g Ryy=%.3g RCN=%.3g delr=%.3f Kor=%.2f Tor=%.2f',nd,J(nd),Kyy(nd),RCN(nd),delr,Kor,Tor));
            lk=3; lw=3; lfp=(lk-1)*lw;
            xdy=[d+(przyr+1+shift)*nDt:nDt:dyf]*dt;
            %xdu=[d1+(przyr+shift)*nDt:nDt:nuf-nDt]*dt;
            if(k==1) subplot(3,1,3); plot(xdy,e,'k'); axis('tight'); xlabel('y-yt');
            else
                subplot(lk,lw,lfp+1); plot(xdy,e,'k'); axis('tight'); xlabel('y-yt');
                subplot(lk,lw,lfp+2); plot(Jx,'k'); axis('tight'); xlabel('Jx');
                subplot(lk,lw,lfp+3); plot(RCNx,'k'); axis('tight'); xlabel('RCNx');
                %subplot(lk,lw,lfp+4); plot(Kyux,'k'); axis('tight'); xlabel('Kyux');
            end
            axx=1;
            %if(abs(delr-0.4)<1.e-4) Ao=A; fio=fi; yo=y;  end
            %input(' ?? ');
        end
    end
    figure(202);
    xd=[1:length(J)]*dt; lk=3; lw=1; lf=lk;
    subplot(lk,1,1); plot(xd,RCN); axis('tight'); xlabel('RCN v.s. opozn.');
    subplot(lk,1,2); plot(xd,J); axis('tight'); xlabel('Wskazzn.MNK v.s. opozn.');
    %subplot(4,1,3); plot(xd,Kyu); axis('tight'); xlabel('Wsp.korelacji Ryu v.s. opozn.');
    subplot(lk,1,lk); plot(xd,Kyy); axis('tight'); xlabel('Wsp.korel. Ryy v.s. opozn.');
    dr=dmin; if(dr<=0) dr=1; end; 
    %nKmin=min(find(Kyu==max(Kyu))); dr=nKmin(1);
    delr=dr*dt; delRC=dRCmax*dt; delNK=dmin*dt;
    % ----- Obl. model dla delR
    d=dr;
    if(1)
        A=Ao; y=yo; yn=yno; Fi=Fio; fi=fio;
    else
        %dyf=d+Ldr; if(dyf>Ld) dyf=Ld; break; end
        %nuf=d1+Ldr-1;
        dyf=d+Ldr; if(dyf>Ld) dyf=Ld; end
        nuf=d1+Ldr-1;
        clear y fi A G yn yp ud;
        %y=Y(d+1:nDt:dyf)'-Y0; % Y(d+1:d+Ldr+1-1)
        %fi(:,1)=Y(d:dyf-1)'-Y0; fi(:,2)=U(d1:nuf)'-U0;
        yn=Y(d+nDt*(shift+1):nDt:dyf)'-y0;
        yp=Y(d+nDt*shift:nDt:dyf-nDt)'-y0; ud=U(d1+nDt*shift:nDt:nuf-nDt)'-u0;
        if(przyr)
            y=diff(yn); yn=y; % Y(d+1:d+Ldr+1-1)
            fi(:,1)=diff(yp); fi(:,2)=diff(ud);
        else
            y=yn; fi(:,1)=yp; fi(:,2)=ud;
        end
        clear yn yp ud;
        lfi=length(y); W=eye(lfi);
        if(wagi==0) fiW=fi'; GW=fi'*fi;
        else
            alf=Ao(1); Kor=Ao(2)/(1-alf); Tor=-DtR/log(alf);
            nWag=round(Tor/DtR);
        end
        nWag=min(nWag,lfi);
        for(i=1:nWag) W(i,i)=(1/(nWag+1-i))^4; end;
        fiW=fi'*W; GW=fiW*fi;
        Ao=inv(GW+R)*fiW*y;
        % -----------------------------------------------
        al=Ao(1); Kor=Ao(2)/(1-al); Tor=-dt/log(al);
    end
end
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
