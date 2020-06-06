% plik KSSident.m
jestDMC=1; idenh=0;
sEhMin=1.e20; sehMin=1.e20; noId=0;
if(jestDMC) m1=1; mf=3; nw1=2; nwf=3; %DMC
else m1=1; mf=2; nw1=5; nwf=6; % PID
end
nf=0; %Lpkth=Ldh;
%for(nw=3:3) %5:5)
for(nw=nw1:nwf)
    if(nw>4) m1=1; mf=2; end% 2:4 wg DMC, 5:.. wg.PID
    for(m=m1:mf)
        nf=nf+1;
        if(nw<5)
            clear txPlik;
            for(i=1:length(typM(m,:)))
                if(typM(m,i)~=' ') txPlik(i)=sprintf('%c',typM(m,i)); end
            end
            txPlik=[txPlik sprintf('W%dRegDMC.csv',nw)];
            plikU=['Pliki/Ur' txPlik]; plikY=['Pliki/Yr' txPlik]; plikYref='Pliki/YrefRegDMC.csv';
            plikDU=['Pliki/DU' txPlik];
            DUr=csvread(plikDU);
        else
            if(m==1) typMod='Skok'; else typMod='Regr'; end %m=3; %=2; % =3, =4;
            txPlik=[Metod(nw-3) typMod 'RegPID.csv'];
            plikU=['Pliki/Ur' txPlik]; plikY=['Pliki/Yr' txPlik]; plikYref='Pliki/YrefRegPID.csv';
            plikDU=['Pliki/DU' txPlik];
            DUr=csvread(plikDU);
        end
        noId=noId+1; 
        %Uid=csvread(plikU); Yid=csvread(plikY); Yref=csvread(plikYref);
        %lU=length(Uid);
        %Ur=Uid(1:nDtR:lU); Yr=Yid(1:nDtR:lU); clear Uid Yid;
        Ur=csvread(plikU); Yr=csvread(plikY); YRef=csvread(plikYref);
        lU=length(Ur);
        %xtu=[1:length(Ur)]*DtR; ntY=[1:length(Yr)]; xty=ntY*DtR;
        %figure(111);
        %plot(xtu,Ur,'b',xty,Yr,'k',xty,YRef(ntY),'r'); axis('tight');
        %DtR=Tos/10; nDt=round(DtR/dt); % okres interwencji regulatora
        %Ld=Ldanych;  Ldr=round(Ld*0.75);
        skok=0; Y0=0; U0=0;  y0=0; u0=0;
        %dMax=2*dr+d1;
        %Ldh=Lpkth;
        hpom=hEmp; 
        Ldhid=length(hpom); Ldh=Ldhid; 
        %while(Ldhid>0.4*Ldh)
            deKonw; %identMNK
            %Ldhid=Ldhid-bezh0-round(1/DtR); 
        %end
    end
end
ident2h;
if(1)
    nDtR=nDt; nDt=1;
    d1=2;
    Ld=length(Y); Ldr=round(Ld*3/4);  %Ldr=round(Ld*0.2);
    dMax=(Ld-Ldr)/3+d1;
    if(dMax>Ldr/2+d1) dMax=floor(Ldr/2)+d1; end
    if(dMax+Ldr>Ld) dMax=Ld-Ldr; end
end
%identMNK;
nDt=nDtR;  Ld=30*tS; %DtR=DtReg; nDt=nDt;



