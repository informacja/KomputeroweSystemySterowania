% RegulPID.m
jaki='p'; txJaki='U'; Metod='iobp'; nwar=0; Model=['Regr';'Skok'];
jaki='d'; txJaki='\DeltaU';
Ur=[]; Yr=[]; YRef=[]; 
for(typR=1:4) %1:4)
    y=[]; X=[]; t=[]; U=[];
    % Nastawy regul.
    if(typR==1 || typR==3)
        del=dels; Kn=Kos; Tn=Tos; %del=1.e13;
        txp='Mod.skok: '; typMod='Skok';
    else
        del=delr; 
        Kn=Kor; Tn=Tor; 
        Kr=1/Kor*Tor/delr*0.2;
        TI=del+Tor; TD=0; fig=6; 
        txp='Mod.regr: '; typMod='Regr';
    end
    for(m=1:4) %m=4;
        nwar=nwar+1; 
        metoda=Metod(m);
        switch(metoda)
            case 'i'
                Kr=1/Kn*Tn/del*0.2;
                TI=del+Tn; TD=0;
            case 'o'
                Kr=1.4*Tn/(Kn*del);
                TI=1.3*del; TD=0.5*del;
            case 'p' % przereg. 20%
                Kr=1.2*Tn/(Kn*del);
                TI=2*del; TD=0.4*del;
            case 'b' % przereg. 20%
                tauR=0.1*Tn; % arbitralne
                Kr=Tn/(Kn*(del+tauR));
                TI=Tn; TD=0; tx=[tx sprintf('T_R=%.2f ',tauR)];
        end
        tx=[txp sprintf('met.%s K_{o}=%.2f T_{o}=%.2f %s=%.2f',metoda,Kn,Tn,txthet,del)];
        %if(jaki=='p') fig=fig+2; end
        % .........................................
        % Symulacja pracy regulatora
        Ur=[]; Yr=[]; YRef=[]; y=[]; U=[]; t=[];; %claear Ur Yr YRef y U t;
        Up=Unom; Upocz=Unom; 
        X=X0; %(1;K1*Up; for(i=2:ntr-invh) X(i)=X(i-1); end 
        %if(invh) X(ntr)=K2*Up; end 
    % Profil wart. zadanej;
        figure(nrSym*4+m); %fig);
        nR=-nDt; 
        Uf=Upocz; dU=0; U0=Upocz; U(1:ntr+2*(nDt-nR))=Upocz;
        U(ntr+2*(nDt-nR):Ld+ntr)=0; 
        y(1:ntr)=Yref(1:ntr); %X(ntr); 
        y2(1:ntr)=X(2); 
        t(1:ntr)=[1:ntr]*dt;
        ep=0; Sp_epp=0; % Stan poczatkowy: 
        yARv1=0; yARv2=0;  Up=Upocz; 
        nkr=0; 
        Yr(nkr+1)=y(ntr); np=1;
        % Przeliczenie ograniczeñ: 
        DUreg=DUhReg*nDt*dt; dUtech=dUhTech*dt; 
        % ............ start zaklocenia ....................
        if(invh<0) v=zeros(1,Ld+1);
        else
            for(n=1:rzad)
                zv=Kv*randn; % rozkl.normalny jest juz standaryzowany
                yARv1=yARv1*alfa+zv; %yARv=yARv1;
                yARv2=yARv2*alfa+yARv1;
                yARv=yARv2;
                v(n)=yARv+srtv;
            end
        end
        % ............. start symulacji dla PID ....................
        if(KSSlab==0 && invh<0) 
            nzFp=0; LdFp=Ld/10; Fpz=Fp; LdFpt=LdFp*rand; 
        end
        for(n=1:Ld)
            ntR=n+ntr; %rzad;
            t(ntR)=ntR*dt;
            if(KSSlab==0) wZakl=0; 
                if(invh<0)
                    wZakl(1)=Tzewn(n); 
                    if(length(vFp)<n)
                        nzFp=nzFp+1;
                        if(LdFpt<=nzFp)
                            LdFpt=LdFp*rand; nzFp=0; Fpz=Fp*(1+0.1*(rand-0.5));
                        end
                        vFp(n)=Fpz;
                    end 
                    wZakl(2)=vFp(n); 
                end
            else
                wZakl=0;
                if(nwar==1)
                    % ------ Zaklocenie wyjsciowe ---------------
                    zv=Kv*randn; % rozkl.normalny jest juz standaryzowany
                    yARv1=yARv1*alfa+zv; %yARv=yARv1;
                    yARv2=yARv2*alfa+yARv1;
                    yARv=yARv2;
                    v(ntR)=yARv+srtv;
                end
            end
            % ------ Koniec symulacji zaklocenia --------
            [X,yn,Umin]=Proces(X,U(n),wZakl);
            Umn(n)=Umin; 
            y(ntR)=yn+v(ntR); % wyjscie
            y2(ntR)=X(2);
            if(nR==nDt) % Regulacja
                Up=U(ntR-1); 
                if(jaki~='p') U0=Up; end
                %[u,du,ep,Sp_epp]=RegPID(y,           Yref,     Kr,TI,TD,Dt, Umin,Umax,DU,   ep,Sp_epp,U0,Up,jaki)
                [Uf,dU,ep,Sp_epp]=RegPID(y(ntR-2:ntR),Yref(ntR),Kr,TI,TD,DtR,Umin,Umax,DUreg,ep,Sp_epp,U0,Up,jaki);
                nkr=nkr+1; 
                DU(nkr)=Uf-Up; Ur(nkr)=Uf; Yr(nkr)=y(ntR); YRef(nkr)=Yref(ntR); vR(nkr)=v(ntR);
                nR=0;
            end
            % ---- Wdrazanie sterowania z uwzgl. ogranicz.przyrostu ---
            % przyrost du jest tu staly (parametr ukl. wykonawczego)
            
            if(abs(dU)>dUtech) dU=sign(dU)*dUtech; end
            if(dU>0 && U(ntR-1)+dU>Uf) dU=Uf-U(ntR-1);
            else if(dU<0 && U(ntR-1)+dU<Uf) dU=Uf-U(ntR-1); end %U(ntR-1)-Uf; end
            %if(dU>0 && U(np)+dU>Uf) dU=Uf-U(np);
            %else if(dU<0 && U(np)+dU<Uf) dU=U(np)-Uf; end
            end
            U(ntR)=U(ntR-1)+dU; if(U(ntR)<Umin) U(ntR)=Umin; end, % Ustawienie sterowania
            if((mod(n,500)==0 || n==Nsym) && KSSlab==0)
                figure(300+nwar);
                subplot(4,1,1); plot(t(1:n),y(1:n),'k',t(1:n),Yref(1:n),'r'); axis('tight'); xlabel('Wyjscie: temperatura w pomieszcz. T(t)');
                subplot(4,1,2); plot(t(1:n),U(1:n),'k',t(1:n),Umn(1:n),'r'); axis('tight'); xlabel('Sterowanie: Temperat. wody goracej Twg(t)');
                vFpMx=max(vFp(1:n)); vFpMn=min(vFp(1:n)); TzewMx=max(Tzewn(1:n)); TzewMn=min(Tzewn(1:n));
                dTzew=TzewMx-TzewMn; dvFp=vFpMx-vFpMn; wFp=dTzew/dvFp;
                subplot(4,1,3); plot(t(1:n),Tzewn(1:n),'k',t(1:n),wFp*(vFp(1:n)-vFpMn),'b'); axis('tight'); xlabel('Zaklocenia: T_{zewn}(t) Fpw(t)');
                subplot(4,1,4); plot(t(1:n),y2(1:n),'k'); axis('tight'); xlabel('Tscian(t) i Twg(t)');
%                 subplot(4,3,11); plot(t,Fwn(1:n),'b'); axis('tight'); xlabel('Fw(t)');
%                 subplot(4,3,12); plot(t,qKa(1:n),'b'); axis('tight'); xlabel('Moc qKa(t)');
                subplot(4,1,1); title(sprintf('Regulacja PID wariant %d metoda=%c',nwar,metoda)); 
            end
            %U(n)=U(np)+dU; % Ustawienie sterowania - TO BYLO NIEPRZYCZYNOWE
            %Up=U(ntR); %Up=U(ntR-1);
            nR=nR+1; np=n;
        end
        figure(nrSym*4+m); %fig);
        %U(nkr+1)=U(nkr); Yref(nkr+1)=Yref(nkr);
        subplot(2,2,typR); hold off;
        if(invh<0) 
            Y0fig=-20; 
        else Y0fig=Y0lin-1; 
        end
        plot(t,U,'b',t,y-Y0fig,'k',t,Yref(1:length(t))-Y0fig,'r'); axis('tight');
        %if(Metod(m)=='o') 
        if(jaki~='p')   
            txPlik=[Metod(m) typMod 'RegPID.csv'];
            plikU=['Pliki/Ur' txPlik]; plikY=['Pliki/Yr' txPlik]; plikYref='Pliki/YrefRegPID.csv';
            plikDU=['Pliki/DU' txPlik]; 
            %lU=length(Ur); csvwrite(plikU,U(1:nDt:lU)); lY=length(y); csvwrite(plikY,y(1:nDt:lY)); lYr=length(Yref); csvwrite(plikYref,Yref(1:nDt:lYr));
            lU=length(Ur); csvwrite(plikU,Ur); lY=length(Yr); csvwrite(plikY,Yr); lYr=length(YRef); csvwrite(plikYref,YRef);
            csvwrite(plikDU,DU);
        end
        ly=length(y); Se=sqrt(sum((y-Yref(1:ly)).^2)/ly); Se=round(Se*1000)/1000; 
        NtSk=[[1:tS] [(tS+2*round(Tos/dt)+ds):tS1]  [(tS1+2*round(Tos/dt)+ds):ly]];
        SE=sqrt(sum((y(NtSk)-Yref(NtSk)).^2)/length(NtSk)); SE=round(SE*1000)/1000;
        txlab=[sprintf('obl.%s %s %st_r=%.2f K_r=%.2f T_I=%.3g T_D=%.2f S_e=%.3g S_{Estab}=%.3g',...
            txJaki,tx,txDel,nDt*dt,Kr,TI,TD,Se,SE)];
        xlabel(txlab);
        subplot(2,2,1); 
        txM='Model z eksperym. '; 
        niel=Nielin+1; txOb=['liniowy  '; 'nielin.1 '; 'nielin.2 '];
        if(nrSym==1) txM=[txM 'czynnego. Obiekt ']; 
        else txM=[txM 'biernego. Obiekt '];
        end
        txM=[txM txOb(niel,:)];
        title(['                                                         ' txM txPar]);
    end
    if(typR==2) jaki='d'; txJaki='\DeltaU'; end
end
subplot(2,2,1); title(['                                                         ' txM txPar]);
