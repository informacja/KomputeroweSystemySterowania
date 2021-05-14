% ===== Robocze ==================
Start=1;

if(Start>0)
    clear all;
    Start=1; 
    if(0) Tm0=15; Tm=20.6+10; To=18.2-3; To0=Tm0; %To=Tm0;
    else  Tm0=0; Tm=1; To0=0; To=0;
    end
    Liz=0.15; % gruboœæ izolacji
    L1=0.3; L2=Liz; L=L1+L2; Ls=L;
    %ParamDyfq;
    if(1) ParamDyf;
    else    
        Spom=121; Hpom=2.8; z=sqrt(Spom); s=4*z; % s=44m;
        Sgs=s*Hpom; % Sgs=44*2.8=123.2 m2 - powierzchnia ogrzewana
        Wzap=14; qwym1m2=Wzap*3600/Sgs; % Wzap [kW] qzap [kW/h/m^2]
        % qwym1m2 moc [kWh/m^2] wymagana na 1m^2 powierzchni scian
        % ......... Parametry cieplne scian ...............
        cws=0.88; ros=1.8e3; lambsc=0.8*3.6; % 0.8 [szac.cw scian [kJ/kg];  ros=1.8e3; szacunk.gêstoœæ [kg/m^3] œciany
        wrcs=ros*cws;
        cwiz=1.460; roiz=40; lambIz=0.035*3.6;
        wrciz=roiz*cwiz;
        cwpow=1.08; ropow=1.23; lambPow=0.025*3.6;
        wrcpow=ropow*cwpow;
        rosz=2500; cwsz=2.5; lambSzk=0.8*3.6;
        wrcsz=rosz*cwsz;
        Diz=lambIz/wrciz; %(roiz*cwiz);
        Dsc=lambsc/wrcs;
        Dszk=lambSzk/wrcsz; %(rosz*cwsz);
        Dpow=lambPow/wrcpow; %(ropow*cwpow);
        D=Dsc;
    end
    % ============== Obliczenia robocze ================
    ifD1nD2=1;
    D1=Dsc; lamb1=lambsc; lamb2=lambIz;
    %wrc1=lamb1/D1; wrc2=lamb2/D2;
    D1=Dsc; D2=Diz; lamb1=lambsc; lamb2=lambIz;
    wrc1=wrcs; wrc2=wrciz;
    if(1)
        Slamb=lamb1*L2+lamb2*L1; wrcSr=(L1*wrc1+L2*wrc2)/Ls;
        D12sr=lamb1*lamb2*Ls/(Slamb*wrcSr);
        wlL1=lamb2*L1/Slamb; wlL2=lamb1*L2/Slamb;
    else
        ifD1nD2=0; D1=Ds; D2=D1; lamb1=lambda; lamb2=lamb1; wrc1=wrcs; wrc2=wrc1;
        D12sr=D1; Slamb=lamb1*Ls; wrcSr=wrc1;
        wlL1=lamb2*L1/Slamb; wlL2=lamb1*L2/Slamb;
    end
    wrc_sr=(wrc1+wrc2);
    Dsr=(lamb1+lamb2)/wrc_sr;
    D1s=lamb1/wrc_sr; D2s=lamb2/wrc_sr;
    % =========================================================
    DT0=Tm0-To0; DT=Tm-To;
    Tg0=(lamb1*Tm0/L1+lamb2*To0/L2)/(lamb1/L1+lamb2/L2);
    qwe=lamb1*DT0/L1;
    clear hTob Txob Txg Txob A B qTx Tx Txob Txgt;
    DL=D1; DR=D1; lambL=lamb1; lambR=lamb1;
    LL=L1; LR=L2; L=Ls;
    if(1)
        ifD1nD2=1;
        D1=Dsc; D2=Diz; lamb1=lambsc; lamb2=lambIz; wrc1=wrcs; wrc2=wrciz;
        % .........................................
        %ifD1nD2=0; D1=Dsc; D2=D1; lamb1=lambsc; lamb2=lamb1; wrc1=wrcs; wrc2=wrc1;
        % .........................................
        wrc_sr=(wrc1+wrc2);
        Dsr=(lamb1+lamb2)/wrc_sr;
        D1s=lamb1/wrc_sr; D2s=lamb2/wrc_sr;
        Slamb=lamb1*L2+lamb2*L1; wrcSr=(L1*wrc1+L2*wrc2)/Ls;
        D12sr=lamb1*lamb2*Ls/(Slamb*wrcSr);
        wlL1=lamb2*L1/Slamb; wlL2=lamb1*L2/Slamb;
    end
    if(ifD1nD2) plik='hTdTg_12.csv'; plikSTg='STg_12.csv'; plikq='hqTL0_12.csv';
    else plik='hTdTg_11.csv'; plikSTg='STg_11.csv'; plikq='hqTL0_11.csv'; end
    hTdTg=csvread(plik);
    Dxx=[0.01 0.0001 5.e-6]; lDx=length(Dxx);
    %tx=[0 100*dt 1000*dt 0.5 1 3 5 6 7 8 9 10 11 12 13 14 15 20 30 50 75 100];
    % .................................................
    DL=D1; DR=D2; LR=L2; L=Ls; mod=1;
    % .................................................
    Nx=1800; Nx1=150; Nx2=300; Nx=3600; Nx=3600000; Nx=4000;
    lf0=2; %lf0=3;
    if(lf0==1) dta=1/3600; dt=dta/10;
    else dta=1/3600; dt=dta; %tx=[0:dt:1/100]; ldt=length(tx);
    end
    tx=[0:dt:100]; lnt=length(tx);
    nxM=lDx; nx0=nxM;
    piL=((pi/L)^2)*DL; wpiL=piL*tx;
    piLR=((pi/LR)^2)*DR; wpiLR=piLR*tx;
    piLL=((pi/LL)^2)*DL; wpiLL=piLL*tx;
    piLL_L=pi*LL/L;
    Int=1;
    %for(Int=1:1) %; %1; %0;
    %for(ndx=nx0:nxM)
    qTLL=zeros(1,lnt); sTg=zeros(1,lnt);
    qTgLR=zeros(1,lnt); qTgLL=zeros(1,lnt);
    hqTg0=zeros(1,lnt); hqTLqL=zeros(1,lnt); % ciep³o od Tg do x=0 i od TLdo 0
    % ....................................
    dxx=Dxx(lDx)*L;
    wxL=dxx/2;
    wXLLdl=pi*wxL/LL; wXLRdl=pi*wxL/LR; XLdl=pi*(LL-wxL)/LL;
    sg=-1; tic;
    for(j=1:Nx)
        j2=j^2;
        sTg=sTg+sin(j*piLL_L)/j*exp(-j2*wpiL); % hTTg T->Tg;
        qTgLR=qTgLR+cos(j*wXLRdl)/j2*exp(-j2*wpiLR); % Tg->q(LL+Dx/2)
        wExp=exp(-j2*wpiLL);
        hqTg0=hqTg0+sg*wExp/j2; hqTLqL=hqTLqL+wExp/j2; % ciep³o od Tg do x=0 i od TLdo 0 
        cosdlexp=wExp*cos(j*wXLLdl)/j2;
        qTgLL=qTgLL+cosdlexp; %exp(-j2*wpiLL); % Tg->q(LL-Dx/2)
        qTLL=qTLL+sg*cosdlexp; %exp(-j2*wpiLL); % T->q(LL-Dx/2)
        %qTLL=qTLL+cos(j*XLdl)/j2*wExp; %exp(-j2*wpiLL); % T->q(LL-Dx/2)
        sg=-sg;
    end
    tim1=toc; fprintf(1,'\nCzas sum=%.2fsek.',tim1);
    dxxpi2=dxx*pi^2;
    Tgxx=LR/L-2/pi*sTg; Tgxx(1)=0; 
    qTgLL=2*LL/(DL*dxxpi2)*qTgLL; qTgLR=2*LR/(DR*dxxpi2)*qTgLR;
    qTLL=2*LL/(DL*dxxpi2)*qTLL;
    DLR=DL*2*lamb1/wrc_sr; DRL=DR*2*lamb2/wrc_sr;
    SqTgTg=(DLR*qTgLL+DRL*qTgLR);
    SqTTg=DLR*qTLL;
    hqTLqL=2*lamb1*LL/(pi^2*DL)*hqTLqL;
    hqTg0=2*lamb1*LL/(pi^2*DL)*hqTg0;
    DbLR=DLR/LL+DRL/LR;
    if(Int)
        Dtx=dta; nDtx=round(Dtx/dt); %(diff(tx));
        lq=length(SqTgTg(1:nDtx:end));
        qTgLL(2:lq)=-DLR*diff(qTgLL(1:nDtx:end))./Dtx; qTgLL(1)=-DLR/LL/dxx;
        qTgLR(2:lq)=-DRL*diff(qTgLR(1:nDtx:end))./Dtx; qTgLR(1)=-DRL/LR/dxx;
        qTLL(2:lq)=-DLR*diff(qTLL(1:nDtx:end))./Dtx; qTLL(1)=-DLR/LL/dxx;
        SqTgTg(2:lq)=-diff(SqTgTg(1:nDtx:end))./Dtx;
        SqTgTg(1)=-(DLR/LL+DRL/LR)/dxx;
        SqTTg(2:lq)=-diff(SqTTg(1:nDtx:end))./Dtx;
        SqTTg(1:4)=-DLR/LL/dxx; %DLR/(LL*dxx)
        hqTg0(2:lq)=lamb1/LL-diff(hqTg0(1:nDtx:end))./Dtx; hqTg0(1)=0;
        hqTLqL(2:lq)=lamb1/LL-diff(hqTLqL(1:nDtx:end))./Dtx; hqTLqL(1)=0;        
    else nDtx=1; Dtx=dt;
    end
    figure(202); 
    subplot(1,2,1); nf1=48; nf2=round(50/dt); 
    plot(tx(1:nf1)*3600,hqTLqL(1:nf1),'k.-'); axis('tight'); xlabel('hqTLqL   t[s]'); 
    subplot(1,2,2); plot(tx(1:nf2),hqTg0(1:nf2),'k'); axis('tight'); xlabel('hqTg0   t[h]');    
    % =================== Teraz temperatury ======================
    Tgx(1)=0; hTdTg(1)=0; Tgx=zeros(1,lnt); dTg1=zeros(1,lnt); dTg2=zeros(1,lnt);
    STg=zeros(1,lnt); ShT0=zeros(1,lnt); hTdTg=zeros(1,lnt); hqTL0=zeros(1,lnt); Tgx=zeros(1,lnt);
    dTgx=[0 diff(Tgxx)];
    tic; opt=1;
    for(n=2:lnt)
        n1=n-1;
        if(1) %opt)
            STg(n)=SqTgTg(2:n)*hTdTg(n-[2:n]+1)';
            ShT0(n)=hqTg0(2:n)*hTdTg(n-[2:n]+1)';
            dTg1(n)=DLR/(LL*dxx)+SqTTg(n);
            dTg2(n)=(DbLR/dxx*Tgx(n1)+STg(n));
        else
            STg(n)=SqTgTg(2:n)*dTgx(n-[2:n]+1)';
            dTg1(n)=DLR/(LL*dxx)+SqTTg(n);
            dTg2(n)=(DbLR/dxx*Tgxx(n1)+STg(n));
        end
        hTdTg(n)=(dTg1(n)-dTg2(n))*Dtx;
        hqTL0(n)=hqTLqL(n)-ShT0(n);         
        Tgx(n)=Tgx(n-1)+hTdTg(n);
    end
    Lnt=n;
    t_obl=toc;
    fprintf(1,'\nCzas=%.2fsek, t.j.: %.0fmin, %.2fsek.',t_obl,floor(t_obl/60),t_obl-floor(t_obl/60)*60);
    csvwrite(plik,hTdTg); csvwrite(plikq,hqTL0); 
else
    if(ifD1nD2) plik='hTdTg_12.csv'; plikSTg='STg_12.csv'; plikq='hqTL0_12.csv';
    else plik='hTdTg_11.csv'; plikSTg='STg_11.csv'; plikq='hqTL0_11.csv'; end
     hTdTg=csvread(plik); Lnt=length(hTdTg);Tgx(Lnt)=0;
     Tgx(1)=0; for(n=2:Lnt) Tgx(n)=Tgx(n-1)+hTdTg(n); end
     Int=1; dt=1/3600; dta=dt; Dtx=dta; nDtx=1; tx=[0:1:Lnt-1]*dt;
end
n1=1; n11=1; nk=100;            wtx=3600;
if(Int) nkt=48/3600/dt; nqt=1; else nkt=12/3600/dt; nqt=10; end
nkk=floor(nkt/nDtx); nkt=floor(nkk*nDtx);
dwa=0; lDX=0;
if(Start>=0)
    figure(lf0+Int);
    kol='b'; kol1='r'; if(nkk<50) kol='b.-'; kol1='r.-'; end
    txDel='\Delta';
    if(Int) txInt=sprintf(' averaged over %st=%.1fs (dt=%.1fs)',txDel,Dtx*3600,dt*3600);
    else txInt=sprintf(' temporary values (dt=%.1fs)',dt*3600); end
    if(0)
        subplot(2,2,1);
        %if(nqt>1&&lf0==1) kol='b'; else kol='b.-'; end
        plot(tx(n1:nDtx:nkt)*wtx,qTgLL(n1:nkk),kol,tx(n1:nDtx:nkt)*wtx,qTgLR(n1:nkk),kol1,tx(n1:nDtx:nkt)*wtx,SqTgTg(n1:nkk),'k');
        %if(nqt>1&&lf0==1) hold on; plot(tx(n1:nDtx*nqt:nkt)*wtx,YQ(n1:nqt:nkk),'b.'); hold off; end
        %text=sprintf('q(t,x)=%sT(t,0-dx)/%sx; %s I_{max}=%d I_{max1}=%d (:)',txDel,txDel,txInt,Nx,Nx1);
        txtS=sprintf('S_{qTbR}(t); %s I_{max}=%d',txInt,Nx);
        texl=sprintf('%sx=%.4f mm',txDel,dxx*1000);
        axis('tight'); xlabel([texl ' t [s]']);  grid on;
        title(txtS);
        subplot(2,2,2);
        nktf=length(tx); nkf=round(length(qTLL)/nDtx);
        plot(tx(n1:nDtx:nktf),qTLL(n1:nkf),'b');
        %if(nqt>1&&lf0==1) hold on; plot(tx(n1:nDtx*nqt:nkt)*wtx,YQ(n1:nqt:nkk),'b.'); hold off; end
        %txtS=sprintf('q(t,x)=%sT(t,0-dx)/%sx; %s I_{max}=%d I_{max1}=%d (:)',txDel,txDel,txInt,Nx,Nx1);
        txtS=sprintf('S_{qTTb}(t); %s I_{max}=%d',txInt,Nx);
        texl=sprintf('%sx=%.4f mm',txDel,dxx*1000);
        axis('tight'); xlabel([texl ' t [s]']);  grid on;
        title(txtS);
        subplot(2,2,3); lw=2; nf=4;
    else subplot(1,2,1); lw=1; nf=2;
    end
    % DLR*qTgLL+DRL*qTgLR
    %plot(tx(n1:nDtx:nkt)*wtx,DRL*qTgLR(n1:nkk),kol,tx(n1:nDtx:nkt)*wtx,DLR*qTgLL(n1:nkk),kol1,tx(n1:nDtx:nkt)*wtx,SqTgTg(n1:nkk),kol1);
    plot(tx(n1:nDtx:nkt)*wtx,SqTgTg(n1:nkk)*dxx,kol);
    texl=sprintf('%sx=%.4f mm',txDel,dxx*1000);
    txtS=sprintf('Heat streams %s %s I_{max}=%d',txInt,texl,Nx);
    texl='q_{TTb}(t)=D_{LR}S_{qTbL}(t)+D_{RL}S_{qTbR}(t)  t[s]';
    axis('tight'); xlabel(texl);  grid on;
    title(txtS);
    subplot(lw,2,nf);
    nktf=length(tx); nkf=round(length(qTLL)/nDtx);
    plot(tx(n1:nDtx:nktf),(DLR/(LL)+SqTTg(n1:nkf)*dxx),'b'), axis('tight');
    xlabel('q_{TL}(t)    t[h]'); grid on;
else opt=1; txDel='\Delta';  
end
figure(27+lf0+opt);
if(ifD1nD2)
    subplot(1,3,1); plot(tx(1:Lnt),Tgx(1:Lnt),'k'); 
    axis('tight'); xlabel('T_{bn}(t)    t[h]');
    subplot(1,3,2); plot(tx(1:Lnt),hTdTg(1:Lnt),'k'); 
    axis('tight'); xlabel(sprintf('h_{TdTb}=%sT_{bL}(t)-%sT_{bR}(t)    t[h]',txDel,txDel));
    
    subplot(1,3,3); plot(tx(1:Lnt),STg(1:Lnt)*Dtx,'k--',tx(1:Lnt),dTg1(1:Lnt)*Dtx,'k',tx(1:Lnt),dTg2(1:Lnt)*Dtx,'k');
    axis('tight'); xlabel(sprintf('Components of h_{TdTb}(t)   t[h]'));
    x1=20+0.01; y1=STg(round(x1/Dtx))*Dtx; x2=10+0.01;  y2=dTg1(round(x2/Dtx))*Dtx; 
    text(x2,y2,sprintf('%sT_{bL}(t),%sT_{bR}(t)',txDel,txDel)); 
    text(x1,y1,'Convolution component'); 
else
    subplot(2,3,1); plot(tx(1:Lnt),Tgx(1:Lnt),'k',tx(1:Lnt),Tgxx(1:Lnt),'r:');
    axis('tight'); xlabel('T_{bn}(t) and T_{ba}(t)    t[h]');
    subplot(2,3,2); plot(tx(1:Lnt),hTdTg(1:Lnt),'k',tx(1:Lnt),dTgx(1:Lnt),'r:');
    axis('tight'); xlabel(sprintf('h_{TdTb}=%sT_{bL}(t)-%sT_{bR}(t)    t[h]',txDel,txDel));
    subplot(2,3,4); plot(tx(1:Lnt),Tgx(1:Lnt)-Tgxx(1:Lnt),'k');
    axis('tight'); xlabel('T_{gn}(t)-T_{ga}(t)    t[h]');
    subplot(2,3,5); plot(tx(1:Lnt),hTdTg(1:Lnt)-dTgx(1:Lnt),'k');
    axis('tight'); xlabel(sprintf('h_{TdTb}(t)-%sT_{ba}(t)    t[h]',txDel));
end
subplot(1,3,3); plot(tx(1:Lnt),STg(1:Lnt)*Dtx,'k--',tx(1:Lnt),dTg1(1:Lnt)*Dtx,'k',tx(1:Lnt),dTg2(1:Lnt)*Dtx,'k');
axis('tight'); xlabel(sprintf('Components of h_{TdTb}(t)   t[h]'));
x1=20+0.01; y1=STg(round(x1/Dtx))*Dtx; x2=10+0.01;  y2=dTg1(round(x2/Dtx))*Dtx; 
text(x2,y2,sprintf('%sT_{bL}(t),%sT_{bR}(t)',txDel,txDel)); 
text(x1,y1,'Convolution component'); 
%end
return;
%end
%end if(mod>1)
