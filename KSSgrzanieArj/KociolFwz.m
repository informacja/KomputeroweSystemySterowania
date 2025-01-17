function [Twg1,Twz1,qkal,Fw,Fwy,brakTwg,rog,roz]=KociolFwz(T,Twg,Twz,qkal);
global wExp CF Cq Kkal_Cq Kkal cicho ro0 cw Fvwz AFgq; 
% gzap>0 charkterystyka qwy(Twg(1:lwg)); 
% gzap<0 obliczenie qwy(Twg); 

% Rownanie wymiennika: tepm.wody Twg i Twz, tem. odbioru T
% Twg-Twz=(Twg-T)*(1-exp(-Kkal/(Fw*cw)));
% Twg-Twg*(1-exp(-Kkal/(Fw*cw)))=Twz-T*(1-exp(-Kkal/(Fw*cw)));
% Twg*(1-1+exp(-Kkal/(Fw*cw)))=Twz-T*(1-exp(-Kkal/(Fw*cw)));
% Twg=(Twz-T*(1-exp(-Kkal/(Fw*cw)))/exp(-Kkal/(Fw*cw));
% qzap=Fw*cw*(Twg-Twz);  
% ............... Dla Fw wymuszonego (staly Fwz): 
% Fw=qFw/(cw*(Twg-Twz))=Fwz*(rog/roz*Twg-Twz)/(Twg-Twz); Fwz - sta�a warto�� przeplywu wody zimnej 
% Fw=Fwz*(rog/roz*Twg-Twz)/(Twg-Twz); 
% roz=ro0/(1+wExp*(Twz-20)); rog=ro0/(1+wExp*(Twg-20));
% Twg=(Twz-T*(1-exp(-Kkal/(Fw*cw)))/exp(-Kkal/(Fw*cw));% ................ Szukamy Twg, kt�re przy ste��j Fwz i Twz daj� zadan� qzap:
% qzap=Fwz*cw*(rog/roz*Twg-Twz);
% St�d 
% 0. roz=ro0/(1+wExp*(Twz-20)); rog=ro0/(1+wExp*(Twg0-20));
% 1. Twg1=(qzap/(Fwz*cw)+Twz)*roz/rog;
% 2. Fw=Fwz*(rog/roz*Twg1-Twz)/(Twg1-Twz); 
% Twg2=T+(Twz-T)/(1-exp(-Kkal/(Fw*cw)));
% 3. Twg2=T+(Twz-T)/(1-exp(-Kkal/(Fw*cw)))
% Z przyr�wnania 1 do 2 mamy
% rog/roz=(1+wExp*(Twz-20))/(1+wExp*(Twg-20))=(1+wExp*(Twg-20)+wExp*(Twz-Twg))/(1+wExp*(Twg-20))
% rog_roz=1-wExp*(Twg-Twz)/(1+wExp*(Twg-20));
% qzap=Fwz*cw*(Twg-Twz)*(1-wExp/(1+wExp*(Twg-20));
brakTwg=0; brak=0; dTmax=5.e-1; 
if(qkal<0) itMax=0; brakM=0; qkal=-qkal; else itMax=400; brakM=1; end
%roz=ro0/(1+wExp*(Twz-20));
iter=0;
while(1)
    iter=iter+1;
    [y,Twg1,Twz1,qwy,Fw,Fwy,rog,roz]=modKotla(T,Twg,Twz,qkal);
    if(iter>=itMax) brakTwg=1; brak=brakM; break; end    
    if(abs(y)<1.e-3) brak=2; break; end
    dTg=1.e-2; 
    Twgn=Twg+dTg;
    y1=modKotla(T,Twgn,Twz,qkal);
    %Twg=Twg-dTg; 
    dy=(y1-y)/dTg; 
    DTg=-y/dy; 
    if(abs(DTg)<1.e-4 && abs(y)<1.e-2) brak=3; break; end
    if(abs(DTg)>dTmax) DTg=sign(DTg)*dTmax; end
    Twg=Twg+DTg; 
end
rog=ro0/(1+wExp*(Twg-20));
qkal=Fvwz*cw*(rog*Twg1-Twz1*roz);
return; 
% ......... Szukamy Twg, ktore da zadana wart. qkal
if(cicho>0) 
    %Twg=AFgq(3)+AFgq(2)*qkal+AFgq(1)*qkal^2; 
    if(qkal<0) TwX=Twg; k0=1;
    else
        Twg=AFgq(2)+AFgq(1)*qkal;
        TwX=[Twg Twg-2:0.001:Twg+2]; k0=1;
    end
    TZ=Twz; EpsTz=1.e-3; LitTz=100;
else k0=2; 
    TwX=[50:0.1:80]; TwX=[Twg-10:0.1:Twg-0.01 Twg Twg+0.01:0.1:Twg+10.01]; 
    %TZ=[35:0.1:49]; TZ=[36:0.1:44]; 
    TZ=[35:0.1:49]; 
    TZ=[30:0.1:70]; TZ=[Twz-10:0.05:Twz-0.05 Twz Twz+0.05:0.05:Twz+10];
    EpsTz=1.e-1; LitTz=1;
end
lTwX=length(TwX); 
% Krok 1. 
% Przy ustalym Twx=TwX(k) szukamy Tz ze wzoru 2: 
lTZ=length(TZ); 
i0=2; iopt=i0; kp=-1;
clear yq; 
for(k=1:lTwX)
    brakTz=0;
    Twx=TwX(k); clear y x DT; 
    rog=ro0/(1+wExp*(Twx-20)); %rog2=rog^2;
    for(i=1:lTZ)
        Tz=TZ(i); 
        for(m=1:LitTz)
            x(i)=Tz;
            mi_sr=lepkw(Twx,Tz,T);
            DT(i)=Twx-Tz;
            if(Tz>Twx-1) 
                y(i)=Tz; brakTz=1; 
                iopt=1; i0=1; k0=1; 
                yq(1)=abs(qkal); yq(2)=yq(1); 
                Fw(k0)=500; 
                break; 
            end
            y(i)=Twx-(Twx-T)*(1-exp(-(Kkal_Cq*mi_sr)/(rog*DT(i)))); Tzob=y(i);
            %Twg-Twz=(Twg-T)*(1-exp(-(Kkal*mi_sr/Cq)/(rog*(Twg-Twz))));
            if(abs(y(i)-x(i))<EpsTz) %=1.e-1)
                i0=i; %Tz=Twx-y(i0);
                if(lTZ==1) Twz=y(i); break; end
            else Tz=y(i); 
            end
        end
    end
    if(brakTz) continue; 
    else
        Twz=y(i0); iopt=i0;
        if(cicho==0)
            figure(300); subplot(2,2,1); 
            plot(x,y,'k',x(i0),y(i0),'r.',x,x,'r'); 
            xlabel(sprintf('Twzo(Tzx) - eq.2: Tzob=%.2f/%.2f Twg=%.2f',x(i0),y(i0),Twx));
            axis('tight'); hold on; 
            subplot(2,2,2); 
            plot(DT,y,'k',DT(i0),y(i0),'r.'); 
            xlabel(sprintf('Twzo(Twg-Tx) - eq.2')); %: Tzob=%.2f/%.2f Twg=%.2f',x(i0),y(i0),Twx));
            axis('tight'); hold on;         
        end
        if(kp<0) kp=k; end
    end
    % Krok 2. 
    % Dla znalezionego Twz=Tz obliczamy qk=yq(k) ze wzoru 1: 
    DTq(k)=Twx-Twz; %DT(i0);
    mi_sr=lepkw(Twx,Twz,T);
    Fw(k)=CF*rog/mi_sr*DTq(k); 
    DTn=(Twg-T)*(1-exp(-Kkal_Cq*CF/(Fw(k))));
    yqk=Cq*rog/mi_sr*DTn^2;
    yq(k)=Cq*rog/mi_sr*(Twx-T)^2*(1-exp(-(Kkal_Cq*mi_sr)/(rog*DTq(k))))^2;
    %wq=[yq(k)/yqk-1, DTn/DTq(k)-1],
    if(qkal<0) k0=k; iopt=i0; break; end
    if(abs(yq(k)/qkal-1)<1.e-2) 
        k0=k; iopt=i0;
        if(cicho>0) break; end
    end
end
qK=yq(k0); Twg=TwX(k0); Twz=y(iopt); 
Fwn=Fw(k0); %CF*(Twg-Twz); 
rog=ro0/(1+wExp*(Twg-20)); roz=ro0/(1+wExp*(Twz-20)); 
if(cicho<=0)
    Q_Tw=(yq(k)-yq(kp))/(TwX(k)-TwX(kp))/1000;
    AFgq=polyfit(yq(kp:k),TwX(kp:k),1); 
    if(cicho==0)
        %yTw=AFgq(3)+AFgq(2)*yq(kp:k)+AFgq(1)*yq(kp:k)^2;
        yTw=AFgq(2)+AFgq(1)*yq(kp:k);
        subplot(2,4,5); plot(TwX(kp:k),yq(kp:k),'k',yTw,yq(kp:k),'m',TwX(k0),yq(k0),'ro');
        axis('tight'); ax=axis; hold on; plot(ax(1:2),[qkal qkal],'r:'); hold off;
        xlabel(sprintf('Moc grzania q(Twg): wzm.Q_Tw=%.2f [kWh/^oC]',Q_Tw));
        subplot(2,4,6); plot(DTq(kp:k),yq(kp:k),'k',DTq(k0),yq(k0),'ro');
        axis('tight'); ax=axis; hold on; plot(ax(1:2),[qkal qkal],'r'); hold off;
        xlabel('       Moc q(Twg-Tz)')
        subplot(2,4,7); plot(TwX(kp:k),TwX(kp:k)-DTq(kp:k),'k',TwX(k0),TwX(k0)-DTq(k0),'ro');
        axis('tight'); xlabel('Twz(Twg)')
        subplot(2,4,8); plot(TwX(kp:k),Fw(kp:k),'k',TwX(k0),Fw(k0),'ro');
        axis('tight'); xlabel('Fw(Twg)')
        subplot(2,2,1); hold off; subplot(2,2,2); hold off;
        subplot(2,4,5); hold off; subplot(2,4,6); hold off;
        subplot(2,4,7); hold off; subplot(2,4,8); hold off;
    end
end