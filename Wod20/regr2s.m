% regr2; 
% =========== Pobieramy dane empiryczne ==========
clear all;  
Ldemp=130; Hammer=0; wspZ=0.02; 
if(Hammer)
    x=[1:Ldemp+1]';
    % alfa obliczymy zak³adaj¹c wart. Nf=Tf/Dt i alfa=exp(-1/Nf); 
    alfa=0.9; 
    Nf=Ldemp/10; alfa=exp(-1/Nf); 
else 
    xmin=0.1; xmax=2.2; 
    dx=(xmax-xmin)/(Ldemp-1); x=[xmin:dx:xmax]'; 
    alfa=0; 
end
if(alfa<1.e-4) Hammer=0; end
[Yemp,Yteor,wspZ,Af]=obiekt(x,-wspZ);
% ......................................................
Kf=length(Af); % Liczba wejsc f.regresji Igo rodz
istot=zeros(1,Kf); % tablica istotnosci wejsc
% ............. Ilustracja danych ......................
if(Hammer)
    figure(1); subplot(1,1,1); plot(x,Yemp,'k.-',x,Yteor,'r'); 
    axis('tight'); xlabel(['\alpha=' sprintf('%.4f',alfa)]); 
    input(' <Ent> '); hold off;
    v0=Yemp(end,1)-Yteor(end,1); 
else
    % ======= Mamy gotowe dane Yemp i liczymy hipotetyczne obserwacje ===
    Ldh=10000; Xhmin=0; Xhmax=xmax*1.2; 
    dx=(Xhmax-Xhmin)/(Ldh-1); xh=[Xhmin:dx:Xhmax]'; 
    [Yh,Yth]=obiekt(xh,wspZ); % Tu wspZ jest odch.stand. zak³
    % ======= Mamy gotowe dane Yemp i x - Rysujemy Y(x) ===
    figure(1); subplot(1,1,1); 
    plot(xh,Yh,'c.',xh,Yth,'b',x,Yemp,'k.'); axis('tight'); 
    input(' <Ent> '); 
end    
% ---------- Budujemy arbitralnie wejscia uogolnione FI ---
tKryt=3.0; % arbitralna wart. graniczna testu Studenta
odrzuc=1; % =0 metoda dolaczania; =1 met.odrzucania
Kd=12; if(Kd>Kf) Kd=Kf; end, % zadana liczba jednomianow (rzad wielom.=(Kd-1)
FId(:,1)=ones(Ldemp,1); istot(1)=1; 
if(Hammer) % Regresja trendu z modelem Hammersteina
    t=x-x(Ldemp+1,1); 
    for(k=2:Kd-1) FId(:,k)=t(2:Ldemp+1).^(k-1); end,
    Ywe=Yemp(1:Ldemp,1); Yewy=Yemp(2:Ldemp+1,1); 
    FId(:,Kd)=Ywe; 
    for(k=2:Kd) istot(k)=1; end, 
    FI=FId; 
else % Regresja statyczna 
    Yewy=Yemp; 
    for(k=2:Kd) FId(:,k)=x.^(k-1); end, % wejscia jednomianowe (arbitralne za³o¿enie)
    if(odrzuc) 
        for(k=2:Kd) istot(k)=1; end, 
    end
    % Standaryzacja wejsc: 
    FIsred=mean(FId); sigFI=[0 std(FId(:,2:Kd))]; 
    FI(:,1)=FId(:,1); Ysred=mean(Yemp); 
    for(n=1:Ldemp) 
        FI(n,2:Kd)=(FId(n,2:Kd)-FIsred(2:Kd))./sigFI(2:Kd); 
    end,
    % ======= Regularyzacja - ocena przydatnosci wejsc =======
    KaraR=[1:Kd-1]*1e-3; % r=1e-3; z regr18
    G=inv([FI(:,2:Kd)'*FI(:,2:Kd)+diag(KaraR)]); 
    Ar=G*FI(:,2:Kd)'*(Yemp-Ysred); % formu³a regresji - wspolcz.modelu
    % ....... Przeliczenie wspolcz. na niestandaryzowane ....
    Am(1,1)=Ysred;
    for(k=2:Kd)
        Am(k,1)=Ar(k-1)/sigFI(k);
        Am(1,1)=Am(1,1)-Am(k,1)*FIsred(k);
    end
    % ------------- Ortogonalizacja -----
[U,S,P]=svd(FI'*FI); % dekomp. singularna: 
%  S macierz diagon. wart.sing. s^2; P- macierz ³aduj¹ca 
U=FI*P; % transformacja ortogonaliz. Karhunena-Loevego: U macierz wejœæ ortogonalnych
s2=0; 
for(i=1:Kd) s2=s2+S(i,i); s(i)=S(i,i); end, % Udzaia³ mocy kolumn ortogonalnych skladowych 

figure(5); plot(s/s2); axis('tight'); 
Aort=inv(U'*U)*U'*(Yemp-Ysred); 
% ------------------------------------------------
    % ........... Wydruk obliczonych wspolcz. ---------
    for(k=1:Kf)
        if(k==1)
            fprintf(1,'\nA0(1)=%-8.3f Am=%-8.3f Aort(1)=%-8.3f .... Af=%-8.3f',...
                         Ysred,Am(1),Ysred,Af(k));
        else
            if(k>Kd) 
               fprintf(1,'\nAr(%2d)=%-8.3f Am=%-8.3f Aort=%-8.3f .... Af=%-8.3f',...
                        k,0,0,0,Af(k));            
               continue; 
            end
            fprintf(1,'\nAr(%2d)=%-8.3f Am=%-8.3f Aort=%-8.3f .... Af=%-8.3f',...
                        k,Ar(k-1),Am(k),Aort(k-1),Af(k));
        end
    end
    input(sprintf('\nModel regularyzowany: KaraR(1)=%.3g <Ent>',KaraR(1))); 
end
clear Am G Ar;  
% ======= Tu zaczynamy regresje krokowa - met. odrzucania 
dalej=1; 
while(dalej)
    clear Fi G Aob Yob E KA Am;  
    kw=1; Fi(:,1)=FI(:,1); 
    for(k=2:Kd) 
        if(istot(k)>0) kw=kw+1; Fi(:,kw)=FI(:,k); end % budujemy Fi tylko z kolumn istotn.
    end
    % ........... Liczymy wspolczynniki modelu .....
    Kw=kw; 
    G=inv(Fi'*Fi); % macierz Gaussa 
    Aob=G*Fi'*(Yewy); % -Ysred); model ze stala A(1)=Ysred% mamy wspolczynniki
    % ..... Policzymy wyjscia modelu dla obserwacji 
    Yob=Fi*Aob; % +Ysred; 
    % ===== Statystyki reszt modelu =========
    E=Yewy-Yob; % reszty (bledy) modelu
    varZ=(E'*E)/(Ldemp-Kw); sigZ=sqrt(varZ); % estym. wariancji Zaklocen Z
    KA=G*varZ;    % macierz kowariancji wspolczynnikow
    % ........... Liczymy wspolczynniki faktyczne Am() ....
    if(Hammer) 
        Am(Kw,1)=Aob(Kw,1); alf=Am(Kw,1); % wspolcz.alfa
        Am(Kw-1,1)=Aob(Kw-1,1)/(1-alf); 
        Am(1,1)=(Aob(1,1)-alf*Am(2,1))/(1-alf); 
        Ysred=Aob(1,1); 
        Afob=Yteor(Ldemp+1); % =Af(1)+Af(2)*(Ldemp+1); 
        % Przeliczenie stalej, bo obiekt liczy czas od 2 do Ldemp+1
        % a regresja ma t=[-Ldemp+2, -Ldemp+3, ..... 0]
        % Dlatego stala Af(1) musi byc przeliczona jw
        kw=Kw; 
        Yf=Fi(:,1:Kd-1)*Am(1:Kd-1,1); 
    else  % ........... Liczymy wspolczynniki faktyczne Am() ....
        kw=1; Am(1,1)=Aob(1); %Ysred; 
        for(k=2:Kd)
           if(istot(k)>0) 
               kw=kw+1; 
               Am(kw,1)=Aob(kw)/sigFI(k); 
               Am(1,1)=Am(1,1)-Am(kw,1)*FIsred(k); 
           end
        end
        Afob=Af(1); 
    end
    % ........... Wydruk obliczonych wspolcz. ---------
    tstmin=1.e40; kw=1; 
    Rgr=tKryt/sqrt(Ldemp-kw+tKryt^2); % graniczny poziom istotn. wsp. korelacji
    Rmax=-1.e20; maxk=1; mink=1000; % indeks wejsc wlaczonych
    tstMale=[]; 
    for(k=1:Kf)
       if(k==1) 
           tSt(k)=abs(Aob(kw))/sqrt(KA(kw,kw)); % liczymt stat. t dla stalej (k=1)
           fprintf(1,'\nAob( 1)=%-8.3f Am=%-8.3f t=%4.2f .... Af=%-8.3f',...
                        Ysred,Am(1),tSt(k),Afob);        
       else
           if(istot(k)<=0) 
               tSt(k)=0; % cz³on nieistotny - kandydat do usuniêcia
               if(k>Kd) Rfe(k)=0;
               else
                   % Liczymy wspolcz.korelacji E z FI(:,k) jeszcze nie w³¹czonymi do modelu  
                   Rfe(k)=(E'*FI(:,k))/(Ldemp*sigZ*std(FI(:,k)));
                   if(Rmax<abs(Rfe(k))) Rmax=abs(Rfe(k)); maxk=k; end
                   if(mink>k) mink=k; end % we.o najnizsz.rzedzie
               end
               fprintf(1,'\nAob(%2d)=%-8.3f Am=%-8.3f R=%4.2f .... Af=%-8.3f',...
                   k,0,0,Rfe(k),Af(k));        
           else 
               kw=kw+1; 
               tSt(k)=abs(Aob(kw))/sqrt(KA(kw,kw)); 
               fprintf(1,'\nAob(%2d)=%-8.3f Am=%-8.3f t=%4.2f .... Af=%-8.3f',...
                            k,Aob(kw),Am(kw),tSt(k),Af(k)); 
               % szukamy min(tSt)
               if(tstmin>tSt(k)) tstmin=tSt(k); kmin=k; end
               if(tKryt>tSt(k)) tstMale=[tstMale, k]; end
           end
       end
    end
    if(Hammer) break; end, % Dla trendu Hammersteina roboczo nie robimy selekcji czlonow istotnych
    kMaxR=max(tstMale); % wejscie o najwyzszym rzedzie  
    if(Kw>1) 
        odrzuc=input(sprintf('\nUsunac <Ent> ; Dodac <1>: Koniec - <-1> ?? ')); 
        if(isempty(odrzuc)) odrzuc=1; 
        else if(odrzuc<0) break; else odrzuc=0; end
        end
    end
    if(odrzuc==0)
        co=input(sprintf('\nDodac we: %d <Ent>(R=%.4f/%.2f) lub <%d> (R=%.4f) lub STOP <-1> ?? ',...
            maxk,Rmax,Rgr,mink,Rfe(mink)));
        if(~isempty(co))
            if(co<0) break; end
            if(co>1 && co <=Kd)
                if(istot(co)>0)
                    fprintf(1,'\n We.%d juz wlaczone',co); continue;
                else maxk=co;
                end
            end
        end
        istot(maxk)=1; fprintf(1,' Wlaczono FI(%d) ',maxk);
    else
        co=input(sprintf('\nUsunac we: %d <Ent>(t=%.4f) lub <%d> (t=%.4f) lub STOP <-1> ?? ',...
            kmin,tSt(kmin),kMaxR,tSt(kMaxR)));
        if(~isempty(co))
            if(co<0) break; end
            if(co>1 && co <=Kd)
                if(istot(co)<=0)
                    fprintf(1,'\n We.%d juz usuniete',co); continue;
                else kmin=co;
                end
            end
        end
        istot(kmin)=0; fprintf(1,' Usunieto FI(%d) ',kmin);
    end
end % dla while(dalej)- koniec regresji krokowej
% .......... odch.stand. bledu estym. i reszt modelu dla obserwacji ......    
LiczbaMalych=0; 
for(n=1:Ldemp) 
   varYm=Fi(n,:)*KA*Fi(n,:)'; varY=varZ+varYm;  
   sgY(n,1)=sqrt(varY);  sgYm(n,1)=sqrt(varYm); 
   if(abs(E(n))<sgY(n)) LiczbaMalych=LiczbaMalych+1; end
end
udzMalych=LiczbaMalych/Ldemp*100; 
if(Hammer==0)
    % .......... odch.stand. bledu estym. i modelu dla dowolnych obserwacji x ......
    Ldv=500; clear fi; 
    dx=(Xhmax-Xhmin)/(Ldv-1); xv=[Xhmin:dx:Xhmax]'; 
    for(n=1:Ldv) 
       % .... wart. wejsc fi dla zadanych xv .....
       fi(1)=1; kw=1; fis(1)=fi(1); 
       for(k=2:Kd) 
           if(istot(k)>0) 
               kw=kw+1; fi(kw)=xv(n)^(k-1); 
               % ........ Standaryz. fi .................           
               fis(kw)=(fi(kw)-FIsred(k))/sigFI(k); 
           end; 
       end
       varYm=fis(1:kw)*KA*fis(1:kw)'; varY=varZ+varYm;  
       sigY(n)=sqrt(varY);  sigYm(n)=sqrt(varYm);  
       Yv(n)=fis(1:kw)*Aob; %+Ysred; 
       Yvo(n)=fi*Am; 
    end
    % ----- Rysujemy wynik ----------
    figure(1); subplot(1,2,1); 
    plot(xh,Yh,'c.',xh,Yth,'b',x,Yemp,'k.',x,Yteor,'k',x,Yob,'r.'), 
    hold on; 
    plot(xv,Yv,'b',xv,Yv+sigYm,'b--',xv,Yv-sigYm,'b--', xv,Yv+sigY,'k--',xv,Yv-sigY,'k--'); axis('tight');
    xlabel(sprintf('Regresja i dane dla Ldemp=%d Kd=%d wspZ=%.3f sigZ=%.3f udzMalych=%.2f %%',Ldemp,Kd,wspZ, sigZ,udzMalych)); 
    hold off;
    subplot(1,2,2); 
    plot(xh,Yh,'c.',xh,Yth,'b',x,Yemp,'k.',x,Yteor,'k',x,Yob,'r.'), 
    hold on; 
    plot(xv,Yv,'b',xv,Yv+sigYm,'b--',xv,Yv-sigYm,'b--', xv,Yv+sigY,'k--',xv,Yv-sigY,'k--'); axis('tight');
    axis('tight'); ax=axis; axis([ax(1) max(x) min([min(Yemp) min(Yteor)]) max([max(Yemp) max(Yteor)])]);
    hold off; 
    %subplot(1,2,2); plot(Yob,E.^2,'k.'); xlabel(sprintf('z^2 v.s. Y_{ob}')); axis('tight'); 
    figure(2); plot(xv,Yv-Yvo); axis('tight'); 
    figure(1);
else % dla Hammersteina - predykcja 
    Ldpr=26; clear fi Ypr; 
    tp=[1:Ldpr]'; xtp=tp+x(Ldemp+1);
    [Ye, Yth]=obiekt(xtp,wspZ,alfa,v0); 
    yn_1=Yemp(Ldemp+1,1); 
    alf=Aob(Kd); 
    % ...... Obliczamy blad trendu dla ost.danej ..........
    vt0=Yemp(Ldemp+1,1)-Am(1,1); 
    for(n=1:Ldpr) 
       % .... wart. wejsc fi dla zadanych xv .....
       fi(1)=1; kw=1; 
       for(k=2:Kd-1) 
           if(istot(k)>0) 
               kw=kw+1; fi(kw)=tp(n)^(k-1); 
           end; 
       end
       fi(Kd)=yn_1; 
       varYm=fi(1:Kd)*KA*fi(1:Kd)'; varY=varZ+varYm;  
       sigYp(n,1)=sqrt(varY);  sigYpm(n,1)=sqrt(varYm);  
       Ypr(n,1)=fi*Aob; 
       yn_1=Ypr(n,1);  
       % teraz wg modelu przeliczonego
       
       v=vt0*alf^n; 
       Yfp(n,1)=fi(1:Kd-1)*Am(1:Kd-1); 
       Yp(n,1)=Yfp(n,1)+v; 
    end
    % ----- Rysujemy wynik ----------
    figure(1); subplot(1,2,1); 
    tx=t(2:Ldemp+1,1);
    plot(x-x(Ldemp+1),Yteor,'c',x-x(Ldemp+1),Yemp,'r.',tx,Yob+sgYm,'r:',tx,Yob-sgYm,'r:',tx,Yob+sgY,'m:',tx,Yob-sgY,'m:',tx,Yob,'b.-',tx,Yf,'k'), 
    hold off; 
    xlabel(sprintf('Regr.Hammerst: Ldemp=%d wspZ=%.0f sigZ=%.0f udzM=%.2f%%',Ldemp,wspZ, sigZ,udzMalych)); 
    %plot(tx,Yob+sigYm(1:Ldemp,1),'b--',tx,Yob-sigYm(1:Ldemp,1),'b--', tx,Yob+sigY(1:Ldemp,1),'k--',tx,Yob-sigY(1:Ldemp,1),'k--'); axis('tight');
    subplot(1,2,2); 
    plot([-1;0;tp],[Yteor(end-1:end,1);Yth(1:end)],'c',0,Yteor(end,1),'c*',[-1;0;tp],[Yemp(end-1:end,1);Ye(1:end)],'b.', 0,Yemp(end,1),'bx');
    hold on; 
    plot(tp(1:end),Ypr-sigYp,'m:',tp(1:end),Ypr+sigYp,'m:',tp(1:end),Ypr-sigYpm,'r:',tp(1:end),Ypr+sigYpm,'r:'); axis('tight');
    plot(tp,Ypr,'b.-',tp,Yp,'g.',tp,Yfp,'k'); 
    hold off;
    xlabel(['Pred: \alpha_o=' sprintf('%.3f/%.3f Lpr=%d',alf,alfa,Ldpr)]); 
    %subplot(1,2,2); 
    %tx=[t(2:Ldemp+1,1); tp];  Yobl=[Yob;Ypr];
    %plot(tx,Yobl,'k',tp,Ye(2:end),'g.',tp,Ypr,'b-.',tx,Yobl+sigYm,'b--',tx,Yobl-sigYm,'b--', tx,Yobl+sigY,'k--',tx,Yobl-sigY,'k--'); axis('tight');
    %hold off; 
end
