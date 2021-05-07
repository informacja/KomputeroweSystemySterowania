% TmRekurencyjne.m
% .................................................................
dtx=dt/5; kt=5;%/5; %/10; 
dtx=dt;
nt=Mdt; tx=[0:nt-1]*dtx; kol='kbmgc'; nrk=2; 
dlx=dlxB; dto=2*dtx;
ParamDyf;
if(0) dlx=dlxB/3; end; %dlx=dlxB; %dlx=0.02; %Ls=L1+L2;
    Lxg=round(L1/dlx)+1; 
    Ldl=round(Ls/dlx)+1; xl=[0:dlx:Ls]; 
    xg=L1; dlx2=dlx*dlx; 
    Tl0(1:Lxg)=(Tm0-Tg0)*(1-xl(1:Lxg)/L1)+Tg0; 
    Tl0(Lxg+1:Ldl)=(Tg0-To0)*(1-(xl(Lxg+1:Ldl)-xl(Lxg))/L2)+To0;     
%else dlx2=dlx*dlx; 
%end
if(1)
    ifD1nD2=1;
    D1=Ds; D2=Diz; lamb1=lambda; lamb2=lambIz; wrc1=wrcs; wrc2=wrciz;
    ifD1nD2=0; D1=Ds; D2=D1; lamb1=lambda; lamb2=lamb1; wrc1=wrcs; wrc2=wrc1;
    wrc_sr=(wrc1+wrc2);
    Dsr=(lamb1+lamb2)/wrc_sr;
    D1s=lamb1/wrc_sr; D2s=lamb2/wrc_sr;
    Slamb=lamb1*L2+lamb2*L1; wrcSr=(L1*wrc1+L2*wrc2)/Ls;
    D12sr=lamb1*lamb2*Ls/(Slamb*wrcSr);
    wlL1=lamb2*L1/Slamb; wlL2=lamb1*L2/Slamb;
    if(0)
        Sxlt1=csvread('Sxlt.csv',20,151); Skxlt1=csvread('Skxlt.csv',20,151);
        Slxlt1=csvread('Slxlt.csv',20,151);
        STg1sr1=csvread('STg1sr.csv'); STg2sr1=csvread('STg2sr.csv');
        SD12sr1=csvread('SD12sr.csv'); SkTgPsr1=csvread('SkTgPsr.csv');
    end
    %dTm=(Tm-Tm0)' 
    %Tf0=(Tf-Tl0)/(Tm-Tm0);
    %b1=-wlL1/L1; b2=-wlL2/L2; bs=-1/Ls; B1=b1-bs; B2=b2-bs; 
    Mm=length(Sxlt(:,1)); 
    for(m=1:Mm)
        %SlB12t(m,Lxg:Ldl)=2/pi*SlB12t(m,Lxg:Ldl)-(Bx2-Bx1)*SEB2(m);
        %hSxl(m,:)=Tf0-2/pi*Sxlt(m,:); %+2/pi*Ssinxt(m,:); 
        %hSxl(m,:)=Tf0-2/pi*(SXl1_2t(m,:)); %-Ssin12xt(m,:)); 
        hSxl(m,:)=Tf0-2/pi*SkLst(m,:)-SlB12t(m,:); %-Ssin12xt(m,:)); 
        hsrSxl(m,:)=Tf0-2/pi*(Skxlt(m,:)-Ssr12xlt(m,:)); 
        h12srSxl(m,:)=Tf0-2/pi*(Sx12lt(m,:)-S12sinxt(m,:)); 
        hsSxl(m,:)=1-xl/Ls-2/pi*Sxlt(m,:); 
        hSlkxl(m,:)=4/pi*(Slxlt(m,:)+0.5*Skxlt(m,:)); 
        hxSlkxl(m,:)=hSxl(m,:)+B1*Ls*2/pi*wlL2*Slxlt(m,:); hxSlkxl(m,Lxg:Ldl)=hxSlkxl(m,Lxg:Ldl)-(B2-B1)*wlL1*Ls*2/pi*Slxlt(m,Lxg:Ldl);
        %[B1*Slxlt(m,1:Lxg-1)+wlL2*Skxlt(m,1:Lxg-1) B2*Slxlt(m,Lxg:Ldl)-wlL2*Skxlt(m,Lxg:Ldl)]; 
        hSlxl(m,:)=-2/pi*[B1*Slxlt(m,1:Lxg-1) B1*Slxlt(m,Lxg)+B2*(Slxlt(m,Lxg:Ldl)-Slxlt(m,Lxg))];
        hnSlxl(m,:)=-[B1*hSlxl(m,1:Lxg-1) B1*hSlxl(m,Lxg)+B2*(hSlxl(m,Lxg:Ldl)-hSlxl(m,Lxg))]; 
    end
    for(m=1:ldt20)  hTX(m,:)=Tf0-[2/pi*(SkLst(m,1:Lxg-1)+SlB12t(m,1:Lxg-1)) Bx2*Ls*(2/pi*(SkLst(m,Lxg:Ldl)+SlB12t(m,Lxg:Ldl)))]; end
    for(m=1:ldt20)  hT(m)=Tf0(Lxg)-2/pi*(SkLxg(m)+SlBxg(m)); wT=(1-hTX(m,Lxg))/(1-hT(m)); hTXs(m,1:Lxg-1)=1-wT*(1-hTX(m,1:Lxg-1)); hTXs(m,Lxg:Ldl)=hTX(m,Lxg:Ldl); end, 
    figure(100); 
    if(0)
        x=[-Ls:0.01:Ls]; Nx=length(x); Skok=zeros(1,Nx); Slin=zeros(1,Nx); 
        sg=-1; 
        for(n=1:20)
            Skok(1,:)=Skok(1,:)+(1-sg)*sin(n*pi*x/Ls)/n*wE;
            Slin(1,:)=Slin(1,:)+sg*sin(n*pi*x/Ls)/n*wE; sg=-sg; 
        end;
        nf=1; figure(119); 
        subplot(2,3,nf); plot(x,2/pi*Skok); axis('tight'); ax=axis; ax(1)=x(1)-0.02; ax(2)=x(end)+0.02; ax(3)=ax(3)-0.02; ax(4)=ax(4)+0.02; axis(ax); xlabel('f_1(x); n_{max}=20'); grid on; 
        subplot(2,3,nf+1); plot(x,2/pi*Slin); axis('tight');ax=axis; ax(1)=x(1)-0.02; ax(2)=x(end)+0.02; ax(3)=ax(3)-0.02; ax(4)=ax(4)+0.02; axis(ax);  xlabel('f_2(x); n_{max}=20'); grid on; 
        subplot(2,3,nf+2); plot(x,2/pi*(Skok+Slin)); axis('tight'); ax=axis; ax(1)=x(1)-0.02; ax(2)=x(end)+0.02; ax(3)=ax(3)-0.02; ax(4)=ax(4)+0.02; axis(ax); grid on; xlabel('f_1(x)+f_2(x); n_{max}=20');
        Skok=zeros(1,Nx); Slin=zeros(1,Nx); 
        for(n=1:1000)
            Skok(1,:)=Skok(1,:)+(1-sg)*sin(n*pi*x/Ls)/n*wE;
            Slin(1,:)=Slin(1,:)+sg*sin(n*pi*x/Ls)/n*wE; sg=-sg; 
        end;
        nf=4;
        subplot(2,3,nf); plot(x,2/pi*Skok); axis('tight'); ax=axis; ax(1)=x(1)-0.02; ax(2)=x(end)+0.02; ax(3)=ax(3)-0.02; ax(4)=ax(4)+0.02; axis(ax); xlabel('f_1(x); n_{max}=20'); grid on; 
        subplot(2,3,nf+1); plot(x,2/pi*Slin); axis('tight');ax=axis; ax(1)=x(1)-0.02; ax(2)=x(end)+0.02; ax(3)=ax(3)-0.02; ax(4)=ax(4)+0.02; axis(ax);  xlabel('f_2(x); n_{max}=20'); grid on; 
        subplot(2,3,nf+2); plot(x,2/pi*(Skok+Slin)); axis('tight'); ax=axis; ax(1)=x(1)-0.02; ax(2)=x(end)+0.02; ax(3)=ax(3)-0.02; ax(4)=ax(4)+0.02; axis(ax); grid on; xlabel('f_1(x)+f_2(x); n_{max}=20');
        
        figure(120);
        x=[-Ls:0.01:Ls]; Nx=length(x); Skok=zeros(20,Nx); Slin=zeros(20,Nx); 
        sg=-1; 
        for(m=1:20) 
            for(n=1:1000) 
                wE=exp(-(n*pi)^2*D1*ta(ndt20(m))); 
                Skok(m,:)=Skok(m,:)+(1-sg)*sin(n*pi*x/Ls)/n*wE; 
                %Slin(m,:)=Slin(m,:)+sg*sin(n*pi*(1-x/Ls))/n*wE; sg=-sg;
                Slin(m,:)=Slin(m,:)+sin(n*pi*(1-x/Ls))/n*wE; sg=-sg; 
            end; 
        end
        nf=1; figure(120); 
        subplot(2,3,nf); plot(x,2/pi*Skok'); axis('tight'); ax=axis; ax(1)=x(1)-0.02; ax(2)=x(end)+0.02; ax(3)=ax(3)-0.02; ax(4)=ax(4)+0.02; axis(ax); xlabel('h_1(x); n_{max}=1000'); grid on; 
        subplot(2,3,nf+1); plot(x,2/pi*Slin'); axis('tight');ax=axis; ax(1)=x(1)-0.02; ax(2)=x(end)+0.02; ax(3)=ax(3)-0.02; ax(4)=ax(4)+0.02; axis(ax);  xlabel('h_2(x); n_{max}=1000'); grid on;
        subplot(2,3,nf+2); plot(x,2/pi*(Skok+Slin)'); axis('tight'); ax=axis; ax(1)=x(1)-0.02; ax(2)=x(end)+0.02; ax(3)=ax(3)-0.02; ax(4)=ax(4)+0.02; axis(ax); grid on; xlabel('h_T=h_1(x)+h_2(x); n_{max}=1000');
        for(m=1:20) hxT(m,:)=1-x/Ls-2/pi*(Skok(m,:)+Slin(m,:)); end; 
        nf=3; subplot(2,2,nf); plot(x(1:end),1-x(1:end)/Ls,'b--',x,hxT(1,:),'k'); axis('tight'); ax=axis; ax(1)=x(1)-0.02; ax(2)=x(end)+0.02; ax(3)=ax(3)-0.02; ax(4)=ax(4)+0.02; axis(ax); grid on; xlabel('Initial and final profile of h_T; n_{max}=1000'); subplot(2,2,nf+1); plot(x(46:end),1-x(46:end)/Ls,'b--',x,hxT(1,:),'r--',x,hxT'); axis('tight'); ax=axis; ax(1)=x(1)-0.02; ax(2)=x(end)+0.02; ax(3)=ax(3)-0.02; ax(4)=ax(4)+0.02; axis(ax); grid on; xlabel('h_T(t,x)=h_1(x)+h_2(x); n_{max}=1000');
        for(m=1:20) hxT(m,:)=x/Ls-2/pi*(Slin(m,:)); end; 
        nf=2; ip=46; subplot(1,3,nf); plot(x(ip:end),x(ip:end)/Ls,'b--',x(ip:end),hxT(1,ip:end),'k',x(ip:end),hxT(:,ip:end)'); axis('tight'); ax=axis; ax(1)=x(ip)-0.02; ax(2)=x(end)+0.02; ax(3)=ax(3)-0.02; ax(4)=ax(4)+0.02; axis(ax); grid on; xlabel('h_{TR}(t,x)'); 
        subplot(1,3,nf+1); plot(x(ip:end),ones(1,46),'b--',x(ip:end),1-1/pi*Skok(:,ip:end)'); axis('tight'); ax=axis; ax(1)=x(ip)-0.02; ax(2)=x(end)+0.02; ax(3)=ax(3)-0.02; ax(4)=ax(4)+0.02; axis(ax); grid on; xlabel('h_{TH}(t,x)');
        % ===== Robocze ==================
        Nx=10000; s=0; s1=zeros(1,length(ta)); 
        for(i=1:Nx) 
            s=s+i*sin(i*pi*LL/L); s1=s1+sin(i*pi*LL/L)/i*exp(-ta*D1*(i*pi/L)^2); 
        end; 
        Tgxx=LR/L-2/pi*s1; figure(2); plot(ta,Tgxx); % Tg
        tx=[100*dt 1000*dt 0.5 1 3 5 6 7 8 9 10 11 12 13 14 15 20 30 50 75 100];
        tx=ta(100:100:end); lnt=length(tx); LR=L2; L=Ls; LL=Ls-LR; k=0; 
        Nx=1000; sTx=zeros(1,lnt); Dxg=0.05; Xg=LL+Dxg; 
        for(j=1:Nx) 
            sTx=sTx+sin(pi*j*Xg/L)/j*exp(-(j*pi/L)^2*D1*tx); 
        end
        Txx=1-Xg/L-2/pi*sTx; 
        Nx=100; 
        z=zeros(1,Nx*Nx); sz=z; sD=z; z1=zeros(1,Nx*Nx); sz1=z1; sD1=z1; ij=0; ji=0; 
        for(j=1:Nx) 
            sinpij=sin(pi*j*Dxg/LR); % dla x-xg=0.05; 
            for(i=1:Nx) 
                x=((i^2*LR^2-j^2*L^2))/((LR*L)^2); 
                if(abs(x)<1.e-60) x=1.e-60; end, 
                if(x<=0) ij=ij+1; z(1,ij)=x; sz(1,ij)=i/j/x*sinpij.*sin(pi*i*LL/L); 
                    sD(1,ij)=i^2*D1*pi^2/L^2; 
                else ji=ji+1; 
                     if(x*D1*tx(lnt)>7.e2) x=7.e2/(D1*tx(lnt)); end 
                     z1(1,ji)=x; sz1(1,ji)=i/j/x*sinpij.*sin(pi*i*LL/L); 
                     sD1(1,ji)=i^2*D1*pi^2/L^2; 
                end
            end; 
        end, 
        y1=zeros(1,lnt); y2=zeros(1,lnt); 
        for(nt=1:lnt) 
            y1(nt)=sum(sz(1,1:ij).*(exp(z(1,1:ij)*D1*tx(nt))-1).*exp(-sD(1,1:ij)*tx(nt)));  
            y2(nt)=sum(sz1(1,1:ji).*(exp(z1(1,1:ji)*D1*tx(nt))-1).*exp(-sD1(1,1:ji)*tx(nt)));  
        end
        dTxn=(4*LR/L^3)/pi^2*y1; dTx2n=(L-Xg)/LR*Tgxx(100:100:end)-Txx;
        %dTx2n=Txx+(L-Xg)/LR*Tgxx(100:100:end)-(L-Xg)/L; 
        figure(1); subplot(2,2,1); hold on; plot(tx,y1,'g'); hold off; axis('tight'); grid on; subplot(2,2,2); hold on; plot(tx,y2,'g'); axis('tight');  hold off; grid on;
        subplot(2,2,3); plot(tx,dTxn,'k',tx,dTx2n,'k'); axis('tight');  grid on; subplot(2,2,4); plot(tx,dTx2n,'g'); axis('tight');  grid on
        
        %xg=LL;
        % Nieco inaczej
        xg=LL;        
        Nx=100; Xg=LL+Dxg; 
        v=zeros(1,Nx*Nx); sU=v; sD=v; v1=zeros(1,Nx*Nx); %sn=v; sn1=v; %sz1=z1; sD1=z1; 
        ij=0; ji=0; 
        for(j=1:Nx) 
            sinpij=sin(pi*j*Dxg/LR)/j; % dla x-x(ip=0.05; 
            for(i=1:Nx) 
                ssin=sinpij*sin(pi*i*LL/L)/i; 
                x=1-(j*L/(i*LR))^2; 
                if(abs(x)<1.e-60) if(sign(x)==0) x=1.e-60; else x=sign(x)*1.e-60; end, end, 
                if(x<=0) ij=ij+1; v(1,ij)=x; sz(1,ij)=ssin/x; 
                    sU(1,ij)=i^2*D1*pi^2/L^2; %sn(1,ij)=ssin;
                else
                    ji=ji+1; sD(1,ji)=i^2*D1*pi^2/L^2;                     
                    %if(x*sD(1,ji)*tx(lnt)>7.e2) x=7.e2/(sD(1,ji)*tx(lnt)); end
                    v1(1,ji)=x; sz1(1,ji)=ssin/x; %sn1(1,ji)=ssin;
                end
            end; 
        end, 
        y1=zeros(1,lnt); y2=zeros(1,lnt); %SU=zeros(lnt,ij); SD=zeros(lnt,ji); 
        sU=sU(1,1:ij); sz=sz(1,1:ij); v=v(1,1:ij); 
        sD=sD(1,1:ji); sz1=sz1(1,1:ji); v1=v1(1,1:ji); 
        for(nt=1:lnt) 
            %SU(nt,1:ij)=sz(1,1:ij).*(exp(v(1,1:ij).*sU(1,1:ij)*tx(nt))-1).*exp(-sU(1,1:ij)*tx(nt));
            %SD(nt,1:ji)=sz1(1,1:ji).*(exp(v1(1,1:ji).*sD(1,1:ji)*tx(nt))-1).*exp(-sD(1,1:ji)*tx(nt));
            %y1(nt)=sum(SU(nt,1:ij));  y2(nt)=sum(SD(nt,1:ji));  
            y1(nt)=sum(sz(1,1:ij).*(exp(v(1,1:ij).*sU(1,1:ij)*tx(nt))-1).*exp(-sU(1,1:ij)*tx(nt)));
            Wx=v1(1,1:ji).*sD(1,1:ji)*tx(nt); nD=find(Wx>7.e2);
            if(length(nD)>0) 
                xd=v1(1,nD); sz1(1,nD)=sz1(1,nD).*xd; xd=7.e2./(sD(1,nD)*tx(nt));
                Wx(nD)=7.e2; sz1(1,nD)=sz1(1,nD)./xd; 
                %if(x*sD(1,ji)*tx(lnt)>7.e2) x=7.e2/(sD(1,ji)*tx(lnt)); end
            end                
            y2(nt)=sum(sz1(1,1:ji).*(exp(Wx)-1).*exp(-sD(1,1:ji)*tx(nt)));  
        end
        y1=y1*(2/pi)^2; y2=y2*(2/pi)^2; 
        % y1 - cz³on wynikaj¹cy ze splotu dTx
        dTx=y1+y2; dTx2=(L-Xg)/LR*Tgxx(100:100:end)-Txx;
        %dTx2n=Txx+(L-Xg)/LR*Tgxx(100:100:end)-(L-Xg)/L; 
        figure(1); subplot(2,2,1); plot(tx,y1); axis('tight'); grid on; subplot(2,2,2); plot(tx,y2); axis('tight');  grid on
        subplot(2,2,3); plot(tx,dTx,'k',tx,dTx2,'r:'); axis('tight');  grid on; subplot(2,2,4); plot(tx,Txx,'k',tx,dTx2+Txx-dTx,'r:'); axis('tight');  grid on

    end
    % ------------- Wariant z Bx1 i Bx2 ----------------------
    if(nowe)
        subplot(1,3,1); plot(xl,Tf0,'b--',xl,2/pi*SkLst+SlB12t); axis('tight');
        subplot(1,3,2); plot(xl,2/pi*SlB12t); axis('tight'); xlabel('Ssin12xt(Ds2/Ls)');
        subplot(1,3,3); plot(xl,(hSxl)); axis('tight'); xlabel('Ssinxt');
    end
    if(nowe==0)
        %subplot(3,3,1); plot(xl,Tf0,'b--',xl,2/pi*(Sxl1_2t-Ssin12xt)); axis('tight');
        subplot(3,3,1); plot(xl,Tf0,'b--',xl,2/pi*SXl1_2t); axis('tight');
        subplot(3,3,2); plot(xl,2/pi*Ssin12xt); axis('tight'); xlabel('Ssin12xt(Ds2/Ls)');
        %subplot(3,3,2); plot(xl,2/pi*Sxlt); axis('tight'); xlabel('Ssin12xt');
        subplot(3,3,3); hold on; plot(xl,(hSxl)); axis('tight'); xlabel('Ssinxt');
        % -------------------------------------
        subplot(3,3,4); plot(xl,Tf0,'b--',xl,2/pi*(Skxlt-Ssr12xlt)); axis('tight');
        subplot(3,3,5); plot(xl,2/pi*Ssinxt); axis('tight'); xlabel('Ssin12xt(Dxsr/Ls)');
        subplot(3,3,6); plot(xl,(hsrSxl)); axis('tight'); xlabel('Ssinxt');
        subplot(3,3,7); plot(xl,Tf0,'b--',xl,2/pi*(Sx12lt-S12sinxt)); axis('tight');
        subplot(3,3,8); plot(xl,2/pi*S12sinxt); axis('tight'); xlabel('Ssin12xt(DsrL12/Ls)');
        subplot(3,3,9); plot(xl,(h12srSxl)); axis('tight'); xlabel('Ssinxt');
    end
    if(1)
        figure(213);
        %subplot(2,4,1); plot(xl,2/pi*[S12xlt(1:Lxg-1) S122xlt(Lxg:Ldl)]');  axis('tight'); ax=axis; axis([[-4 Ldl+4]*dlxB ax(3:4)]); xlabel('Sxlt');
        %subplot(2,4,2); plot(xl,2/pi*Ssr12xlt'); axis('tight'); ax=axis; axis([[-4 Ldl+4]*dlxB ax(3:4)]); xlabel('Slxlt');
        subplot(2,4,3); plot(xl,2/pi*Ssinxt'); axis('tight'); ax=axis; axis([[-4 Ldl+4]*dlxB ax(3:4)]); xlabel('Skxlt');
        subplot(2,4,1); plot(xl,2/pi*S12xlt');  axis('tight'); ax=axis; axis([[-4 Ldl+4]*dlxB ax(3:4)]); xlabel('Sxlt');
        subplot(2,4,2); plot(xl,2/pi*Skxlt');  axis('tight'); ax=axis; axis([[-4 Ldl+4]*dlxB ax(3:4)]); xlabel('Sxlt');
        subplot(2,4,4); plot(xl,2/pi*Ssin12xt');   ax=axis; axis([[-4 Ldl+4]*dlxB ax(3:4)]); xlabel('Slxlt+0.5*Skxlt');
        
        subplot(2,3,4); plot(xl,hsSxl');  axis('tight'); ax=axis; axis([[-4 Ldl+4]*dlxB ax(3:4)]); xlabel('Sxlt');
        subplot(2,3,5); plot(xl,hxSlkxl',xl,Tf0,'b--'); axis('tight'); ax=axis; axis([[-4 Ldl+4]*dlx ax(3:4)]); xlabel('hSlxlt');
        subplot(2,3,6); plot(xl,hSxl'); axis('tight'); ax=axis; axis([[-4 Ldl+4]*dlxB ax(3:4)]); xlabel('hSxl');
    %end
    else
        figure(213);
        subplot(2,4,1); plot(2/pi*Sxlt');  axis('tight'); ax=axis; axis([-4 Ldl+4 ax(3:4)]); xlabel('Sxlt');
        subplot(2,4,2); plot(2/pi*Slxlt'); axis('tight'); ax=axis; axis([-4 Ldl+4 ax(3:4)]); xlabel('Slxlt');
        subplot(2,4,3); plot(2/pi*Skxlt'); axis('tight'); ax=axis; axis([-4 Ldl+4 ax(3:4)]); xlabel('Skxlt');
        subplot(2,4,4); plot(-hSlkxl');   ax=axis; axis([-4 Ldl+4 ax(3:4)]); xlabel('Slxlt+0.5*Skxlt');
        subplot(2,3,4); plot(hsSxl');  axis('tight'); ax=axis; axis([-4 Ldl+4 ax(3:4)]); xlabel('Sxlt');
        subplot(2,3,5); plot(xl,hxSlkxl',xl,Tf0,'b--'); axis('tight'); ax=axis; axis([[-4 Ldl+4]*dlx ax(3:4)]); xlabel('hSlxlt');
        subplot(2,3,6); plot(hSxl'); axis('tight'); ax=axis; axis([-4 Ldl+4 ax(3:4)]); xlabel('hSxl');
        %subplot(2,4,4); plot(-2/pi*(Skxlt'+Sxlt'));   ax=axis; axis([-2 Ldl+2 ax(3:4)]); xlabel('Skxlt+Sxlt');
    end
else D12sr=D1; Slamb=lamb1*Ls; wrcSr=wrc1;
end
alf1=exp(-D1*dto/dlx2); alf2=exp(-D2*dto/dlx2); alfs=exp(-Dsr*dto/dlx2);
alf1_1p2=(1-alf1)/2;  alf2_1p2=(1-alf2)/2;
alfs1_1p2=(1-alfs)*(D1s/Dsr); alfs2_1p2=(1-alfs)*(D2s/Dsr);
clear A B Tx Txg; 
A=zeros(Ldl-2,Ldl-2); 
if(1)
    a11=alf1; a12=alf1_1p2;  a22=alf2; a21=alf2_1p2; 
    as=alfs; as1=alfs1_1p2;  as2=alfs2_1p2; 
else
    a11=1-D1*dto/dlx2; a12=D1*dto/dlx2/2;  a22=1-D2*dto/dlx2; a21=D2*dto/dlx2/2;
    as=1-Dsr*dto/dlx2; as1=D1s*dto/dlx2;  as2=D2s*dto/dlx2; 
end
A(1,1)=a11; A(1,2)=a12; 
if(0) 
    A(Ldl-2,Ldl-2)=a11; A(Ldl-2,Ldl-3)=a12; 
    for(m=2:Ldl-3) A(m,m)=a11; A(m,m-1)=a12; A(m,m+1)=a12; end
else
    A(Ldl-2,Ldl-2)=a22; A(Ldl-2,Ldl-3)=a21; 
    for(m=2:Lxg-2) A(m,m)=a11; A(m,m-1)=a12; A(m,m+1)=a12; end
    m=Lxg-1; A(m,m)=as; A(m,m-1)=as1; if(m<Ldl-2) A(m,m+1)=as2; end
    for(m=Lxg:Ldl-3) A(m,m)=a22; A(m,m-1)=a21; A(m,m+1)=a21; end
end
%C=inv(A); %B=exp(A*dtx); C=(B.^nt); 
%Tx=Tl0(2:Ldl-1)'; 
%Tx=C*Tl0(2:Ldl-1)'+(eye(Ldl-2)-C)*[Tm*D1/2/dlx2;zeros(Ldl-4,1);To0*D1/2/dlx2];  Tx(Lxg-1)/Txgt(nt,4)-1, [Tx(Lxg-1) Txgt(nt,4)],
tic; Tx=Tl0(2:Ldl-1)'; 
clear hTob Txob;
Txg=zeros(1,nt); qTx=zeros(1,nt); %Txgt=zeros(nt,6);
Txob=zeros(ldt20,Ldl); hTob=zeros(ldt20,Ldl); 
B=[Tm*a12;zeros(Ldl-4,1);To0*a21]; 
%Tx=zeros(Ldl-2,1); B=[(Tm-Tm0)*a12;zeros(Ldl-4,1);(To0-To0)*a21];
mx=1; 
for(n=1:nt) 
    Tx=A*Tx+B; 
    Txg(n)=Tx(Lxg-1); 
    qTx(n)=lamb1*(Tm-Tx(1))/dlx;
    if(abs(tx(n)-ta(ndt20(mx)))<=dtx) 
        Txob(mx,:)=[Tm Tx' To0]; hTob(mx,:)=(Txob(mx,:)-Tl0)/(Tm-Tm0); 
        if(mx<ldt20) mx=mx+1; end, 
    end
    Txgt(n,1:6)=[Tx([1 Lxg-3:Lxg+1])]';
end 
figure(103); 
subplot(1,3,1); plot(xl,Txob); axis('tight'); grid; 
subplot(1,3,2); plot(hTob'); axis('tight'); grid; 
subplot(1,3,3); plot(hSxl'); axis('tight'); ax=axis; ax(3)=0; axis(ax); grid;
q0=lamb1*(Tm0-Tg0)/L1; 
hqTx=(qTx-q0)/(Tm-Tm0); 
Dt=5/3600; ndtx=round(Dt/dtx); ndt=round(Dt/dt); 
nxf=length(hqTx); N=round(nt/ndtx); m=0; for(n=1:ndtx:nxf-ndtx-1) m=m+1; hQTx(m)=mean(hqTx(n:n+ndtx-1)); end
nf=length(hqTm); N=round(nt/ndt); m=0; for(n=1:ndt:nf-ndt-1) m=m+1; hQTm(m)=mean(hqTm(n:n+ndt-1)); end
figure(211);  subplot(1,2,1); %hold on;
n1=1; tf=20/3600; n2=round(tf/dtx); m2=round(tf/dt); 
lqTx=length(hQTx); lqTm=length(hQTm); nq2=round(tf/Dt); %mq2=round(tf/Dt); 
plot(tx(1:n2)*3600,hqTx(1:n2),'k',tx(1:ndtx:n2+1)*3600,hQTx(1:nq2+1),'k--',...
    t(1:m2)*3600,hqTm(1:m2),'r',t(1:ndt:m2+1)*3600,hQTm(1:nq2+1),'r--',...
    [tx(1) tx(n2)]*3600,[hqTx(end) hqTx(end)],'k',t(1:m2)*3600,hqTm(1:m2),'r',...
    [tx(1) tx(n2)]*3600,[hqTm(end) hqTm(end)],'r:',...
    tx(1)*3600,hqTx(1),'k.',t(1)*3600,hqTm(1),'r.',t(1)-dt,0,'w'); grid on;
xlabel('hqTx(t) k rekur.; hqTm(t) r: analit. t[sek]'); axis('tight'); hold off;
subplot(1,2,2); dnt=round(dt/dtx);
%plot(tx(n1:dnt:n2)*3600,hqTx(n1:dnt:n2)./hqTm(n1:m2)-1,'k',[tx(1) tx(n2)]*3600,[hqTx(end)./hqTm(end)-1 hqTx(end)./hqTm(end)-1],'k:'); xlabel('Blad wzgledny hqTx/hqTm-1'); axis('tight'); grid on;
plot([1:nq2+1]*Dt*3600,hQTx(1:nq2+1)-hQTm(1:nq2+1),'k',[1 nq2+1]*Dt*3600,[hQTx(end)-hQTm(end) hQTx(end)-hQTm(end)],'k:'); xlabel('Blad wzgledny hqTx/hqTm-1'); axis('tight'); grid on;
nxt=length(hTg); 
%Tx=Tx+Tl0(2:Ldl-1)'; Txg=Txg+Tl0(Lxg-1); Tx0gt=[Tl0([1 Lxg-3:Lxg+1])]'; for(k=1:6) Txgt(:,k)=Txgt(:,k)+Tx0gt(k); end
%Tx(Lxg-1)/Txgt(nt,4)-1, [Tx(Lxg-1) Txgt(nt,4)],
fprintf(1,'\nTmRek: Czas=%.3g sek.\n',toc);
Tgt=(Tm-Tm0)*hTg+Tg0; 
Tgt2=Tgf-(Tm-Tm0)*2/pi*(lamb1*Ls*SD2'-lamb1*L2*SkTgP')/(lamb1*L2+lamb2*L1); 
if(1)
    Tw1=(Tm-To0)*[lamb2*L1/(lamb1*L2+lamb2*L1) lamb1*L2/(lamb1*L2+lamb2*L1) lamb1*Ls/(lamb1*L2+lamb2*L1)];
    [Tw1;  Tm-Tgf Tgf-To0 Tm-To0]
  
    Tx00=(Tm-To0)*(-2/pi*(lamb1*L2*Skxl'))/(lamb1*L2+lamb2*L1); 
    Tx01=(Tm-Tgf)*(1-xlb'/L1)-(Tm-To0)*(2/pi*(lamb2*L1*Sxl-lamb1*L2*Skxl)')/(lamb1*L2+lamb2*L1);
    figure(223); subplot(1,2,1); 
    plot([0:length(Skxl)-1]*dlxB, -2/pi*(Sxl-Skxl),'k', [0:length(Sxl)-1]*dlxB, -2/pi*Skxl(1:end),'k--',...
        [0:length(Sxl)-1]*dlxB, -2/pi*Sxl(1:end),'m'); 
    subplot(1,2,2); 
    plot([0:length(Tx00)-1]*dlxB, Tx00,'k--',[0:length(Tx01)-1]*dlxB, Tx01,'k'),
end
Tm-Tm0
figure(212);  subplot(1,2,1); %hold on;
plot(tx(1:nt),Txg(1:nt),kol(1),t(1:nxt),Tgt(1:nxt),'r:',[tx(1) tx(nt)],[Txg(nt) Txg(nt)],'k:',[t(1) t(nxt)],[Tgt(nxt) Tgt(nxt)],'r:'); grid on;
xlabel('Txg(t) k rekur.; Tgt(t) r analit. t[godz]'); axis('tight'); hold off;
n1=1; subplot(1,2,2); plot(tx(n1:nxt),Txg(n1:nxt)./Tgt(n1:nxt)-1,'k',[tx(1) tx(nxt)],[Txg(nxt)./Tgt(nxt)-1 Txg(nxt)./Tgt(nxt)-1],'k:'); xlabel('Blad wzgl.Txg(t)/Tgt(t)-1. t[godz]'); axis('tight'); grid on;
if(ifD1nD2==0)
    Tgt=Txgo(:,2)'; Tgt(1)=Tgt(2); dTg=diff(Tgt); dTg=[dTg dTg(end)]; 
    if(ifD1nD2) nrTg2=3; TxgL=Tf(Lxg-1)-(Tm-Tm0)*2/pi*SL1; else TxgL=Txgo(:,1)'; nrTg2=4; end
    TxgL(1)=TxgL(2); Txgo(1,nrTg2)=Txgo(2,nrTg2); 
    % Txgo(:,1) dla xg-dl; Txgo(:,2)=Tgt dla xg; Txgo(:,4) dla xg+dl; 
    dT=dt/dlxB*(D1*TxgL-Dsr*Tgt+D2*Txgo(:,nrTg2)')/dlxB;
    figure(215); 
    subplot(2,2,1); plot(t,dTg,'k',t,dT,'r:'); xlabel('dTg i d^2T/dl^2'); axis('tight'); subplot(2,2,3); plot(t,dTg-dT,'k'); xlabel('dTg-d^2T/dl^2');  axis('tight');
    dxT=dto/dlx*(D1s*Txgt(:,3)'-Dsr*Txg+D2s*Txgt(:,5)')/dlx;
    dxTgt=diff(Txg); dxTgt=[dxTgt dxTgt(end)]; 
    subplot(2,2,2); plot(tx,dxTgt,'k',tx,dxT,'r:'); xlabel('dxTg(k) i d^2Tx/dl^2 r:'); axis('tight'); subplot(2,2,4); plot(tx,dxTgt-dxT,'k'); xlabel('dxTg-d^2Tx/dl^2');  axis('tight');
end
%subplot(1,2,2); hold on; 
%plot(tx(1:nt),Txg(1:nt)-Tgt(1:nt),'r:'); 
%xlabel('Txg(rekur.)-Tgt(analit.)'); axis('tight'); 
Tlx=[Tm Tx' To0];
%return; 
Tl=Tl0; Tl(1)=Tm; Tl(Ldl)=To0; 
nf1=102; figure(nf1); subplot(1,2,1);
plot(xl(1:Ldl),Tl0,'k',xl(1:Ldl),Tf,'r'), hold on; %xp,'m',xl(1:Ldl),Txp2,'b',xl(1:Ldl),Txp4,'k');
tic;
for(n=1:Mdt) %Mdt)
    %dT_dt=D1/2*(Tl(1:Lxg-2)-2*Tl(2:Lxg-1)+Tl(3:Lxg));
    %Tl(2:Lxg-1)=Tl(2:Lxg-1)+dT_dt*dt_dlx2;
    Tl(2:Lxg-1)=Tl(2:Lxg-1)*alf1+alf1_1p2*(Tl(1:Lxg-2)+Tl(3:Lxg));
    % Tl(2,n+1)=Tl(2,n)+D1/2*(Tl(1,n)-2*Tl(2,n)+Tl(3,n));
    % Tl(2,n+1)=Tl(2,n)-Tl(2,n)*D1+D1/2*(Tl(1,n)+Tl(3,n));
    % ................................................
    %  Dla sekwencji k, k+1 plastrów mamy
    % Tl(k+1,n)=Tl(k,n)+dl*dT(k,n)/dl
    % Tl(k+1,n+1)=Tl(k,n+1)+dl*(T(k+1,n+1)-T(k,n+1))/dl)
    % Tl(k,n+1)=Tl(k,n)*alf+(1-alf)/2*(Tl(k-1,n)+Tl(k+1,n));
    % Tl(k+1,n+1)=Tl(k+1,n)*alf+(1-alf)/2*(Tl(k,n)+Tl(k+2,n));
    % Tl(k,n+1)=(Tl(k-1,n)+dl*dT(k-1,n)/dl)*alf+(1-alf)/2*(Tl(k-1,n)+Tl(k+1,n));
    % Tl(k,n+1)=(Tl(k-1,n)*(1+alf)/2+(1-alf)/2*Tl(k+1,n))+dl*dT(k-1,n)/dl*alf;
    % ................................................
    %dT_dt=-Dsr*Tl(Lxg)+Dsr*(D1s/Dsr*Tl(Lxg-1)+D2s/Dsr*Tl(Lxg+1));
    %Tl(Lxg)=Tl(Lxg)+dT_dt*dt_dlx2;
    Tl(Lxg)=Tl(Lxg)*alfs+alfs1_1p2*Tl(Lxg-1)+alfs2_1p2*Tl(Lxg+1);
    % ................................................
    %dT_dt=D2/2*(Tl(Lxg:Ldl-2)-2*Tl(Lxg+1:Ldl-1)+Tl(Lxg+2:Ldl));
    %Tl(Lxg+1:Ldl-1)=Tl(Lxg+1:Ldl-1)+dT_dt*dt_dlx2;
    Tl(Lxg+1:Ldl-1)=Tl(Lxg+1:Ldl-1)*alf2+alf2_1p2*(Tl(Lxg:Ldl-2)+Tl(Lxg+2:Ldl)); ;
    % ........................................
    qT0(n)=lamb1*(Tl(1)-Tl(2))/dlx;
    Txgt(n,:)=[Tl(2) Tl(Lxg-2:Lxg+2)];
    if(t(n)<100 && mod(n,round(Mdt/50))==0 || t(n)>=100 && mod(n,round(Mdt/5))==0)
        figure(nf1);
        subplot(1,2,1);
        plot(xl(1:Ldl),Tl,'k'), hold on; %xp,'m',xl(1:Ldl),Txp2,'b',xl(1:Ldl),Txp4,'k');
        axis('tight'); %xlabel('Profil temperatury dla ta=');
        subplot(1,2,2); nlp=2;
        plot(t(np:n),Txgt(np:n,nlp),'b:',t(np:n),Txgt(np:n,nlp+1),'b-.',t(np:n),Txgt(np:n,nlp+2),'k',t(np:n),Txgt(np:n,nlp+3),'r-.',t(np:n),Txgt(np:n,nlp+4),'r:'), %xp,'m',xl(np:Ldl),Txp2,'b',xl(np:Ldl),Txp4,'k');
        axis('tight'); hold on; xlabel(sprintf('n=%d t=%.2f',n,t(n)));
        aa=1;
    end
    if(n==nt1) Tx1=Tl;
    else if(n==ntp) Txp=Tl;
        else if(n==ntp2) Txp2=Tl;
            else if(n==ntp4) Txp4=Tl;
                end
            end
        end
    end
end
fprintf(1,'\nCzas rek.=%.3g sek.\n',toc);
figure(nf1); subplot(1,2,1); hold off; subplot(1,2,2); hold off;
figure(101);
subplot(1,3,1);
plot(xl(1:Ldl),Tx1,'k',xl(1:Ldl),Tl,'k',xl(1:Ldl),Txp,'k',xl(1:Ldl),Txp2,'k',xl(1:Ldl),Tl0,'r',xl(1:Ldl),Tf,'r',xl(1:Ldl),Txp4,'k');
axis('tight'); xlabel('Profil temperatury dla t='); hold off;
subplot(1,3,2);
if(normTxg)
    for(k=-2:0) Txgt(:,4+k)=(Txgt(:,4+k)-(Tm0-Tg0)*(1-xl(Lxg+k)/L1)-Tg0)/DTm; end
    for(k=1:2)  Txgt(:,4+k)=(Txgt(:,4+k)-(Tg0-To0)*(1-(xl(Lxg+k)-xl(Lxg))/L2)-To0)/DTm; end
end
plot(t,Txgt(:,2),'k:',t,Txgt(:,3),'k--',t,Txgt(:,4),'k',t,Txgt(:,5),'k-.',t,Txgt(:,6),'k:');
axis('tight'); ; xlabel('Odpowiedzi skokowe Tl dla xg-2dl: xg i xg+2dl'); hold off;
subplot(1,3,3);
plot(t,qT0,'k'); axis('tight'); xlabel('qT0'); hold off;
csvwrite('Txgt.csv',Txgt(:,[1 3:5]),Mdt,4);
