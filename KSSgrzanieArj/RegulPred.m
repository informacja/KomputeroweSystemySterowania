%RegulPred
global jestF111 typReg;
jestF111=0;
nPrdYref=5*Nu; nPrdYref=1; %Nopt;
wagaR=[1 49 121 225]; lW=length(wagaR);
typM=['skokowy';'regres.';'pomiar ']; % model: s skokowy,r regr, p z pomiaru
Du=Umax-Unom;
for(nw=2:2) %1:lW)
    for(m=1:3) %:-1:1)
        % ======= Synteza regulat. predykcyjnego DMC ======
        Nu=5; wR=wagaR(nw); 
        Ko=Kos; To=Tos; del=round(dels/dt);
        if(typM(m,1)=='r') Ko=Kor; To=Tor; del=round(delr/dt); end
        if(typM(m,1)=='p') Ko=Kob; U(1)=Unom; Del=0; To=Tos+dels; end
        [WR,h,t,Ihmax,rDu,DelR,Nopt]=SyntDMC(typM(m,1),U(1),Du,Ko,To,del,vh(nSkok:end),DtR,Nu,wR);
        if(typM(m,1)=='p') hpom=h; end
        if(jestF111) 
            typReg=sprintf('Reg.DMC Model %s W%d delR=%.3f',typM(m,:),nw,DelR*DtR);
        end
        %[WR,h,t,Ihmax,rDu,DelR,Nopt]=SyntDMC(metoda,U0,Du,Ko,To,do,DtR,Nu,wR)
        if(nw==1)
            figure(11);
            subplot(3,1,m); plot(t,h,'b.-'); axis('tight'); txDelta='\Delta'; txthet='\theta';
            xlabel(sprintf('Odp.skok: %st_R=%.2f i %sU=%.2f %s_h=%.2f Ihmax=%d',txDelta,DtR,txDelta,Du,txthet,DelR*DtR,Ihmax));
        end
        % teraz liczymy odpowiedz
        clear y X t U DU t y U;
        Upocz=Unom; %Yref(1)/Kob;
        Uf=Upocz; dU=0; U0=Upocz; U(1:ntr)=Upocz;
        U(ntr+1:Ld+ntr)=0; nR=0; %-nDt; 
        DU=zeros(1,round(Ld/nDt)); %Ur=zeros(1,round(Ld/nDt));
        t(1:ntr)=[1:ntr]*dt; 
        Yref(Ld+ntr+1:Ld+ntr+nPrdYref*nDt)=Yref(Ld+ntr);
        % Stan poczatkowy: 
        Up=Upocz; X=X0;
        %X(1)=K1*Up; for(i=2:ntr-invh) X(i)=X(i-1); end 
        %if(invh) X(ntr)=K2*Up; end
        y(1:ntr)=Yref(1); 
        y2(1:ntr)=X(2);         
        nkr=1; nr=0;
        clear Ur Yr YRef;
        Yr(nkr)=y(ntr); Ur(nkr)=Up; DU(nkr)=0; Ufp=Up;
        for(n=1:Ld)
            ntR=n+ntr; %rzad;
            t(ntR)=ntR*dt;
            % ------ Zaklocenie wyjsciowe v(n) jak dla PID ---------------
            wZakl=0; if(KSSlab==0 && invh<0) wZakl(1)=Tzewn(n); wZakl(2)=Fp; end
            [X,yn,Umin]=Proces(X,U(n),wZakl);
            y(ntR)=yn+v(ntR); % wyjscie
            y2(ntR)=X(2); 
            if(nR==nDt) % Regulacja
                nr=nr+ntr; nkr=nkr+1; 
                Up=U(ntR-1); DU(nkr)=0; Ur(nkr)=Up;
                n0=n-round(Ihmax*DtR); if(n0<1) n0=1; end
                Y0=h(Ihmax)*U(n0);
                dUf=RegDMC(Yref(n:nDt:n+nPrdYref*nDt)',WR,h,y(ntR),Nopt,DU(1:nkr),Ur(1:nkr),DelR,Umax,Umin, DUreg, nkr);
                %DU=RegDMC(Yref,WR,h,yn,Nopt,dU, Y0, DelR, Up, Umax, Umin, DUreg,n)
                dU=dUf(1,1);
                Uf=Up+dU;  DU(nkr)=dU;
                Ur(nkr)=Uf; Yr(nkr)=y(ntR); YRef(nkr)=Yref(ntR);
                nR=0;
                if(mod(nkr,30)==0 && jestF111)
                    lU=length(Ur);
                    xtu=[1:length(Ur)]*DtR; ntY=[1:length(Yr)]; xty=ntY*DtR;
                    figure(111);
                    plot(xtu,Ur,'b',xty,Yr,'k',xty,YRef(ntY),'r'); axis('tight');
                end
                if(~jestF111 && (mod(n,500)==0 || n==Ld) && KSSlab==0)
                    figure(400+m);
                    subplot(4,1,1); plot(t(1:n),y(1:n),'k',t(1:n),Yref(1:n),'r'); axis('tight'); xlabel('Wyjscie: temperatura w pomieszcz. T(t)');
                    subplot(4,1,2); plot(t(1:n),U(1:n),'k'); axis('tight'); xlabel('Sterowanie: Temperat. wody goracej Twg(t)');
                    subplot(4,1,3); plot(t(1:n),Tzewn(1:n),'b'); axis('tight'); xlabel('Zaklocenie: T_{zewn}(t)');
                    subplot(4,1,4); plot(t(1:n),y2(1:n),'k'); axis('tight'); xlabel('Tscian(t) i Twg(t)');                    
                end
            end
            subplot(4,1,1); title(sprintf('Regulacja PREDYKCYJNA waga=%.3f model %s',wagaR(nw),typM(m,:))); 
            % ---- Wdrazanie sterowania z uwzgl. ogranicz.przyrostu ---
            % przyrost du jest tu staly (parametr ukl. wykonawczego)
            %dUtech=DU_tech*1000; 
            if(abs(dU)>dUtech) 
                dU=sign(dU)*dUtech; 
            end
            if(dU>0 && U(ntR-1)+dU>Uf) dU=Uf-U(ntR-1);
            else if(dU<0 && U(ntR-1)+dU<Uf) dU=Uf-U(ntR-1); end %U(ntR-1)-Uf; end                
%            else if(dU<0 && U(ntR-1)+dU<Uf) dU=U(ntR-1)-Uf; end
            end
            U(ntR)=U(ntR-1)+dU; % Ustawienie sterowania
            nR=nR+1;
        end
        %U(nkr+1)=U(nkr); Yref(nkr+1)=Yref(nkr);
        if(nw>1)
            clear txPlik;
            for(i=1:length(typM(m,:))) 
                if(typM(m,i)~=' ') txPlik(i)=sprintf('%c',typM(m,i)); end
            end
            txPlik=[txPlik sprintf('W%dRegDMC.csv',nw)];
            plikU=['Pliki/Ur' txPlik]; plikY=['Pliki/Yr' txPlik]; plikYref='Pliki/YrefRegDMC.csv';
            plikDU=['Pliki/DU' txPlik]; 
            %lU=length(Ur); csvwrite(plikU,U(1:nDt:lU)); lY=length(y); csvwrite(plikY,y(1:nDt:lY)); lYr=length(Yref); csvwrite(plikYref,Yref(1:nDt:lYr));
            lU=length(Ur); csvwrite(plikU,Ur); lY=length(Yr); csvwrite(plikY,Yr); 
            csvwrite(plikDU,DU); lYr=length(YRef); csvwrite(plikYref,YRef);
        end
        figure(nrSym*10);
        nf=(nw-1)*3+m;
        subplot(lW,3,nf); Lt=length(t);
        plot(t(1:Lt),U(1:Lt),'b',t,y(1:Lt)-Y0fig,'k',t,Yref(1:Lt)-Y0fig,'r'); axis('tight');
        ly=length(y); Se=sqrt(sum((y-Yref(1:ly)).^2)/ly); Se=round(Se*1000)/1000; 
        SE=sqrt(sum((y(NtSk)-Yref(NtSk)).^2)/length(NtSk)); SE=round(SE*1000)/1000;
        txlab=['Mod.' sprintf('%s rDu=%.1f %st_r=%.2f Nu=%d Nopt=%d S_e=%.3g S_{Estab}=%.3g',...
             typM(m,:),rDu,txDel,DtR,Nu,Nopt,Se,SE)];
        xlabel(txlab);
        if(jestF111) 
            figure(111); hold off; 
        end        
    end
end
figure(nrSym*10); subplot(lW,3,2); title([' Regulator predykcyjny. Obiekt ' txOb(niel,:) sprintf('rzad %d',rzad)]); %                                     Obiekt rzedu ' txPar]);