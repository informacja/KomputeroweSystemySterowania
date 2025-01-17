% SymulGrzanie.m
%global ndrg ndrz; 
Tn(1)=T; Tsn(1)=Ts; Twg(1)=Twgn; qKa(1)=qKaln; Twzn(1)=Twz;
%Tsym=48; Tsym=12; 
Nsym=round(Tsym/dt); t=[]; t(1)=0; 
if(Fwconst==0)
    NDel=ceil(Lrg/dlmin); pozFwg(1:Ndel)=[1:Ndel]*dlg; pozFwz(1:Ndel)=[1:Ndel]*dlg; 
end
%ntSk=round(8/dt); DU=5; 
Tzew=Tzewn(1); TwgX=ones(1,nTwg1)*Twg(1); 
for(n=2:Nsym)
    t(n)=(n-1)*dt;
    if(Fwconst==0) dlg=Fwn(n-1)/sr*dt/rog; end
    nd=n-ndrg; if(nd<1) nd=1; end, % opozn. wody goracej
    [Tp(n) Tn(n) Tsn(n) Twzn(n) TwgX Fwn(n) qKa(n) rog roz]=Ogrzew(Tn(n-1),Tsn(n-1),Twzn(nd),TwgX,Twg(nd),Tzew,-qKaln,Fp); 
    if(Fwconst==0) 
        dlz=Fwn(n)/sr*dt/roz; pozFwz(1:Ndel)=pozFwz(1:Ndel)+dlz; 
        dlg=Fwn(n)/sr*dt/rog; pozFwg(1:Ndel)=pozFwg(1:Ndel)+dlg; 
    end
    nd=n-ndrz; if(nd<1) nd=1; end % opoznienie wody zimnej
    nd=n-1; 
    Twg(n)=Twg(nd); %if(n<round(1/dt)) Tzew=Tzewn(n); end
    if(n==ntSk) %round(8/dt))% piec musi uzupelnic warto�� qKa(n)
        Twg(n)=Twg(nd)-DU; %Tzew=Tzewn(n);
    end
    Tzew=Tzewn(n);
    %dlg=Fwn(n)/sr*dt/rog; pozFwg(1:Ndel)=pozFwg(1:Ndel)+dlg; 
    if(mod(n,300)==0 || n==Nsym)
        figure(300);
        subplot(4,1,1); plot(t,Tn(1:n),'k'); axis('tight'); xlabel('Temperatura w pomieszcz. T(t)');
        subplot(4,1,2); plot(t,Tsn(1:n),'r'); axis('tight'); xlabel('Temperat. scian Ts(t)');
        subplot(4,1,3); plot(t,Tzewn(1:n),'b'); axis('tight'); xlabel('T_{zewn}(t)');        
        subplot(4,3,10); plot(t,Twzn(1:n),'b',t,Twg(1:n),'r'); axis('tight'); xlabel('Twz(t) i Twg(t)');
        subplot(4,3,11); plot(t,Fwn(1:n),'b'); axis('tight'); xlabel('Fw(t)');
        subplot(4,3,12); plot(t,qKa(1:n),'b'); axis('tight'); xlabel('Moc qKa(t)');
    end
end
