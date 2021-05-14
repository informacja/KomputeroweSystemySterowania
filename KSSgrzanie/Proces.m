function [X, Y, Umin]=Proces(X,U,w)
global rzad rz1 Tob K1 K2 dt invh WspXn C1 Cf S1 Sf;
global lX nTwg nTwg1; 
%K1=Kob^(1/rzad); Ko=K1^(rzad-1); K2=K1*(2/(rzad-1));
if(invh<0) 
   %[Tp(n) Tn(n) Tsn(n) Twzn(n) Twgwy Fwn(n) qKa(n) rog roz]=Ogrzew(Tn(n-1),Tsn(n-1),Twg(nd),Twzn(nd),Tzew,-qKaln,Fp); 
    [Y X(1) X(2) X(3) X(4:lX)]=Ogrzew(X(1),X(2),X(3),X(4:lX),U, w(1),-1,w(2)); 
    % X(lX+[1:nTwg1]) - Temperat.wody zimnej w rurach
    X(lX+[1:nTwg])=X(lX+[2:nTwg1]); X(lX+nTwg1)=X(3); 
    Umin=max([X(1),w(1),4,X(lX+1)]); %min(X(1),X(3)); Umin=T; %min(T,Twz); 
    return;
end
Umin=0; 
if(WspXn<0)
    W=U*K1+w;
    for(i=1:rzad)
        if(i<rzad || invh==0) dX=-X(i)/Tob+W/Tob;
        else dX=-X(i)/Tob+(K2*U+w)/Tob;
        end
        X(i)=X(i)+dX*dt; % Calkowanie metod¹ Eulera
        W=X(i);
    end
    if(invh) Y=X(i-1)-X(i);  %ostatni cz³on generuje nieminim.
    else Y=X(i);
    end
else
    %K1=F/(4*S1); C1=c1/S1; u nalezy do (0,2) Umax=1.8
    %K2=F/Sf; Cf=cf/Sf
    % u=U0+U/4; u0=0.5; 4*u0=2.0; 
    [F1,F2]=przepl(U);
    W=F1/S1+w; 
    rzad1=rz1-1;
    for(i=1:rzad1)
        Fwy=C1*sqrt(X(i)-WspXn*X(i+1)); 
        dX(i)=W-Fwy;
        W=Fwy;
    end
    for(i=1:rzad1)
        X(i)=X(i)+dX(i)*dt; % Calkowanie metod¹ Eulera
        if(X(i)<1.e-2) X(i)=1.e-2; end
    end
    i=rz1;
    Fwy=C1*sqrt(X(i)); 
    X(i)=X(i)+(W-Fwy)*dt; if(X(i)<1.e-2) X(i)=1.e-2; end
    Y=C1*S1*sqrt(X(i)); 
    if(invh) 
        i=rzad;
        %W=F2; %K2=F/(4*Sf); Cf=cf/Sf
        Fwy=Cf*sqrt(X(i)); 
        X(i)=X(i)+(F2/Sf-Fwy)*dt; if(X(i)<1.e-2) X(i)=1.e-2; end
        Y=Y+Cf*Sf*sqrt(X(i)); 
    end
end
