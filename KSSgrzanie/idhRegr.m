% idhRegr.m
function [Kor,Tor,dr,delR,Ndrf]=idhRegr(Y,U,y0,u0,dt,nDt,DtR,d1,dMax,Ldr,Ld,nskok,jestK)
global przyr wagi PokazId; 
nd=0; 
tSkok=nskok*dt; 
clear RCN J Je Jx Kyu Kyy;
Jdmin=1.e20; Jdmax=0; RCNmax=-1; dRCmax=-1; 
Yd=Y;
if(jestK) Kor=jestK; R=0; kf=1; 
else  kf=2;    R=[2 0; 0 1]; R=[1 0; 0 1]*0.;
end
% d1=1; Tak powinno byÊ zawsze
npl=find(Y(d1:nDt:end)*sign(Y(end))>0); 
if(~isempty(npl)) d01=(npl(1)-1)*nDt+1; else d01=1; end
% yd(1)=Y(1)=yh1(nskok-nDt)  yd(2)=Y(nDt+1)=yh1(nskok-nDt+nDt)=yh1(nskok)=0; 
% yd(3)=Y(1+2*nDt)=yh1(nskok-nDt+2*nDt)=yh1(nskok+nDt) mÛze byÊ wiÍksze od 0   
% Y(1+2*nDt) odpowiada zewn. wartoúci nDt po chwili skoku (Del=1); jesli npl(1)=1 
%       tzn. d01=2*nDt, tzn.del=nDt;, a wiÍc mamy obiek úciúle w≥aúciwy bez opoznienia
% d1 to 1.sza prÛbka ciπg≥a w okresie dyzkretyzacji nDt
nsk_1=nskok-1; kmin=1; dmin=d01; delmin=d01-2*nDt;
d=d01-nDt; Koniec=0; ndmin=2; 
if(d1+Ldr>Ld-15*nDt) Ldr=Ld-d1-15*nDt; end 
if(d+Ldr>Ld-10*nDt) Ldr=Ld-d-10*nDt; end 
while(d<dMax)
    d=d+nDt;
%for(d=d01:nDt:dMax)
    %Koniec=0;
    if(d+Ldr>Ld) 
        if(Koniec==3) break; end
        Koniec=3; d=dmin-nDt; nd=ndmin-1; kkmin=kmin; 
        continue; 
    else dyf=d+Ldr; 
        if(Koniec==0) 
            ndf=nd+1; 
        end
    end
    nuf=d1+Ldr-1; nuf=d1+Ldr+nDt; % bylo +1
    nd=nd+1; %nd=d-d1+1;    
    clear y Fi Fid fi A G yn yp ud Dud;
    %y0=Y(d); u0=U(d1);
    del=d-2*nDt; dy1=del+nDt;%(del=d-nDt+nDt; bo zaczynamy od -nDt
    if(1)
        % Y(nsk_1) jest wartoﬂci· Y dla kroku poprzedzaj·cego skok (Y(nsk_1)=0 dla pojedyncz. skoku 
        yd=[Y(d:nDt:dyf)]'; %-y0;
        %yp=[Y(nsk_1) Y(d-nDt:nDt:dyf-nDt)]'-y0; 
        yp=[Y(nsk_1) Y(d:nDt:dyf-nDt)]'; %-y0; 
        ud=U(d1:nDt:nuf-nDt)'; %-u0;
    else
        % Y(nsk_1) jest wartoﬂci· Y dla kroku poprzedzaj·cego skok (Y(nsk_1)=0 dla pojedyncz. skoku 
        yd=[Y(nsk_1) Y(d:nDt:dyf)]'; %-y0;
        %yp=[Y(nsk_1) Y(d-nDt:nDt:dyf-nDt)]'-y0; 
        yp=[Y(nsk_1) Y(nsk_1) Y(d:nDt:dyf-nDt)]'; %-y0; 
        ud=U(d1:nDt:nuf)'; %-nDt)'; %-u0;
    end        
    if(przyr)
        y=diff(yd); yd=y; % Y(d+1:d+Ldr+1-1)
        Fid(:,1)=diff(yp); Fid(:,2)=diff(ud); 
    else y=yd; Fid(:,1)=yp; Fid(:,2)=ud; 
    end
    clear yp yn fi Fi A; 
    fi(:,1:kf)=Fid(:,1:kf); 
    if(jestK) 
        fi(:,1)=fi(:,1)-Kor*ud; y=y-Kor*ud; A(2,1)=Kor; 
    end
    Fi=fi; yn=y;
    kt=1; deltR=floor(del/nDt);
    km=round(deltR/4); if(km==0) km=1; end
    km=2;
    clear y fi xyt;
    lyn=length(yn(:,1));
    obsGap=0; 
    k1=1; kkonc=km; 
    obsGap=0; 
    if(Koniec && obsGap) k1=kmin; kkonc=kmin; end
    for(k=k1:kkonc)
        clear lykk y fi xyt fiu;
        kk=k;
        if(kk>1 && obsGap)
            lykk=length(yn(kk:end,1));
            y=[yn(1,1)';yn(kk:end,1)]; fi=[Fi(1,:);Fi(kk:end,:)]; 
            if(jestK) fiu=Kor*[Fid(1,2); Fid(kk:end,2)]; end
            xyt=[d [d+(kk-1)*nDt:nDt:d+(kk-2+lykk)*nDt]]*dt;
        else
            y=yn(1:end,1); fi=Fi(1:end,:); 
            if(jestK) fiu=Fid(1:end,2); end
            xyt=[d:nDt:d+(lyn-1)*nDt]*dt;
        end            
        lfi=length(fi(:,1)); W=eye(lfi);
        if(0) %k>1)
            fi(:,1)=yt(1:lfi);
        end
        if(wagi==0) fiW=fi'; GW=fi'*fi;
        else
            if(k>1)
                alf=A(1); Tor=-DtR/log(alf);
                if(jestK==0) Kor=A(2,1)/(1-alf); else A(2,1)=Kor*(1-alf); end
                nWag=round(Tor/DtR);
            else nWag=deltR;
            end
            nWag=min(nWag,lfi);
            W(1,1)=0; wykl=1; %wykl=2; wykl=4; 
            for(i=2:nWag) W(i,i)=(1/(nWag+1-i))^wykl; end;
            fiW=fi'*W; GW=fiW*fi;
        end
        S=svd(GW); RCNx(k)=S(end)/S(1);
        A(1:kf,1)=inv(GW+R)*fiW*y; %if(A(1)>1) continue; end
        alf=A(1,1); 
        if(jestK) A(2,1)=Kor*(1-alf); else Kor=A(2,1)/(1-alf); end
        clear yt; yt=fi*A(1:kf,1);
        %for(i=1:lfi) yt(i,1)=yp*A(1)+fi(i,2)*A(2); yp=yt(i,1); end
        ysr=mean(y); ypsr=mean(fi(:,1));
        usr=mean(Fid(kk:end,2)); 
        e=y-yt;
        Jx(k)=(e'*e)/length(yt);
        %Kyux(k)=((y-ysr)'*(fi(:,2)-usr))/(length(y));
    end
    Jxk(nd)=Jx(k); 
    shift=k-1;
    if(dt>=0.99999) delr=del*DtR; else delr=del*dt; end
    alf=A(1,1); Tr=-DtR/log(alf); if(Tr<1.e-6) continue; end
    Tor=Tr;
    if(jestK==0) Kor=A(2,1)/(1-alf); else A(2,1)=Kor*(1-alf); end
    %clear yn Fi;
    %J(nd)=(e'*W*e)/length(yt);
    %if(delr>2*Tor) break; end
    RCN(nd)=S(end)/S(1);
    if(RCN(nd)>RCNmax) RCNmax=RCN(nd); dRCmax=nd; end
    if(przyr)
        sgy=std(y); sgu=std(Fid(kk:end,2)); sgyp=std(fi(:,1));
    else sgy=1; sgu=1; sgyp=1;
    end
    %Kyu(nd)=((y-ysr)'*(fi(:,2)-usr))/(length(y)*sgy*sgu);
    Kyy(nd)=((y-ysr)'*(fi(:,1)-ypsr))/(length(y)*sgy*sgyp);
    % ======= Obliczamy wyjﬂcie ====================
    lfi=length(fi);
    clear ytm ytf;
    % obliczamy model przekszta?cony: y_n-Ku_n-d=alf*(y_n-1-Ku_n-d-1)
    ytm(1,1)=0; ynp=ytm(1,1); dFid=0; %fi(1,1); 
    for(i=2:lfi)
        if(jestK) ynp=ytm(i-1,1)+dFid; 
        else ynp=ytm(i-1,1)+fi(i-1,2)*A(2,1); %ytm(i-1,1)=ynp;
        end
        ytm(i,1)=ynp*A(1,1);
        if(i<lfi && jestK) dFid=fiu(i-1,1)-fiu(i,1); end
        %if(jestK==0) ytm(i,1)=ytm(i,1)+fi(i,2)*A(2,1); end
        %ynp=ytm(i,1); 
    end
    % .................................................
    % obliczamy model oryginalny: y_n=alf*y_n-1-(1-alf)*Ku_n-d-1
    LFi=length(Fid(:,1)); ynp=Fid(1,1);
    for(i=1:LFi)
        ytf(i)=ynp*A(1,1)+Fid(i,2)*A(2,1);
        ynp=ytf(i);
    end
    if(Koniec) break; end
    % Obliczamy wskanik
    Em=ytm-y;
    Jm(nd)=(Em'*Em)/length(ytm-kf);
    E=yd(2:end)-ytf(2:end)'; sigE=sqrt(E'*E/length(E-kf-1));
    JE(nd)=sigE;
    if(1) J(nd)=JE(nd); else J(nd)=Jm(nd); end
    if(J(nd)>Jdmax) Jdmax=J(nd); delmax=del; end; %Jdmin=1.e20;
    if(J(nd)<Jdmin)
        delmin=del; Jdmin=J(nd); dmin=d; 
        Ao=A; fio=fi; Fio=Fi; yno=yn; yo=y; kmin=kk; ndmin=nd;
    else
        if(Jdmin>0 && delr>2*Tor)
            Koniec=Koniec*10+1; d=dmin-nDt; nd=ndmin-1; kkmin=kmin;
            continue; % break;
        end
        if(0) %nd-delmin>10 && delmax<delmin)
            if(Jdmax/Jdmin>1.50) 
                Koniec=Koniec*10+2; d=dmin-nDt; nd=ndmin-1; kkmin=kmin;
                continue; %break; 
            end
        end
    end
    if(1) 
        plotfig200; end; %mod(nd,nDt)==0)
    if(abs(yd(2)/Y(end))>0.5 && Jdmin>0 && Koniec==0) 
        Koniec=Koniec+4; d=dmin-nDt; nd=ndmin-1; kkmin=kmin; continue; 
    end
    if(Koniec) break; end
end
if(Koniec)
    dr=delmin; d=dmin; kk=kmin; Ndrf=Ldr+d; 
    delRC=dRCmax*dt; delNK=delmin*dt;
    % ----- Obl. model dla delR
    del=dr; delR=round(dr/nDt)*DtR;
    A=Ao; y=yo; yn=yno; Fi=Fio; fi=fio; kk=kmin; nd=ndmin; 
    alf=Ao(1,1); Kor=Ao(2,1)/(1-alf); Tor=-DtR/log(alf);    
    fprintf(1,'\nKoniec=%d nd=%d dr=%d delR=%.2f',Koniec,nd,dr,delR); 
    % ======= Obliczamy wyjﬂcie ====================
    lfi=length(fi);
    clear ytm ytf;
    % obliczamy model przekszta?cony: y_n-Ku_n-d=alf*(y_n-1-Ku_n-d-1)
    ytm(1,1)=0; ynp=ytm(1,1); dFid=0; %fi(1,1);
    for(i=2:lfi)
        if(jestK) ynp=ytm(i-1,1)+dFid;
        else ynp=ytm(i-1,1)+fi(i-1,2)*A(2,1); ytm(i-1,1)=ynp;
        end
        ytm(i,1)=ynp*A(1,1);
        if(i<lfi) dFid=Fid(kk+i-2,2)-Fid(kk+i-1,2); end
        %if(jestK==0) ytm(i,1)=ytm(i,1)+fi(i,2)*A(2,1); end
        %ynp=ytm(i,1);
    end
    % .................................................
    % obliczamy model oryginalny: y_n=alf*y_n-1-(1-alf)*Ku_n-d-1
    LFi=length(Fid(:,1)); ynp=Fid(1,1);
    for(i=1:LFi)
        ytf(i)=ynp*A(1,1)+Fid(i,2)*A(2,1);
        ynp=ytf(i);
    end
    if(1) plotfig200; end; %mod(nd,nDt)==0)
end
if(0)
    figure(2);
    xd=[1:length(J)]*dt; lk=3; lw=1; lf=lk;
    subplot(lk,1,1); plot([1:length(RCN)],RCN); axis('tight'); xlabel('RCN v.s. opozn.');
    subplot(lk,1,2); plot(xd,J); axis('tight'); xlabel('Wskazzn.MNK v.s. opozn.');
    %subplot(4,1,3); plot(xd,Kyu); axis('tight'); xlabel('Wsp.korelacji Ryu v.s. opozn.');
    subplot(lk,1,lk); plot([1:length(Kyy)],Kyy); axis('tight'); xlabel('Wsp.korel. Ryy v.s. opozn.');
end

