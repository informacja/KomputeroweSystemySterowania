% plik SyntDMC
function [WR,h,t,Ihmax,rDu,DelR,Nopt]=SyntDMC(metoda,U0,Du,Ko,To,do,vh,DtR,Nu,wR)  
global Unom X0 Y0lin WspXn rzad Tob Kob K1 K2 dt invh;
global hPom hEmp;
Y0=Ko*U0; h(1)=0; 
DelR=ceil(do*dt/DtR+1); 
ntR=round(DtR/dt);
Nopt=round(4*To/DtR); 
Ldmin=(Nopt+DelR)*ntR; 
Ldh=10*Ldmin; Ihmax=Ldh; 
U(1:Ldh)=Unom+Du; nR=0; i=1; hp=0; h(1)=0; yp=Y0;
if(metoda~='p') % na podstawie modelu
    rDu=wR*Nopt/Nu*Ko^2;
    for(n=do+1:Ldh)
        y=model(Ko,To,do,U,n,yp);
        nR=nR+1; yp=y;
        if(nR==ntR && y/Ko>0.001)
            i=i+1; nR=1;
            h(i)=(y-Y0)/Du;
            t(i)=(n-do-1)*dt;
            if(n>Ldmin && abs((h(i)-hp)/Ko)<0.001)
                Ihmax=i;
                break;
            end
            hp=h(i); %(hp+h(i))/i;
        end
    end
    aa=i;
else % z pomiaru
    DelR=1; Nopt=0;
    wgPom=1;
    if(wgPom==0)
        Y0=Y0+Y0lin;
        if(WspXn<0)
            X(1)=K1*U0; for(i=2:rzad-invh) X(i)=X(i-1); end
            if(invh) X(rzad)=K2*U0; end        
        else X=X0;
        end
        vn=vh(1); i=0;
    end
    for(n=1:Ldh)
        if(wgPom==0)
            [X y]=Proces(X,U(n),0);
            if(n<length(vh)) vn=vh(n+1); end
            y=y+vn;
        end
        nR=nR+1; Lhp=length(hPom); 
        if(nR==ntR)
            i=i+1; nR=1;
            if(wgPom) if(n<=Lhp) h(i)=hPom(n); else h(i)=h(i-1); end
            else h(i)=(y-Y0)/Du; 
            end
            t(i)=n*dt;
            if(n>Ldmin)
                if((h(i)-hp)/Ko<0.01 && Nopt==0) Nopt=i;
                else if((h(i)-hp)/Ko<0.001) Ihmax=i; break; end
                end
            end
            hp=h(i);
        end
    end
    if(wgPom) hEmp=h; end
end
rDu=wR*Nopt/Nu*Ko^2; 
% Synteza macierzy H
H=zeros(Nopt,Nu); 
for(nu=1:Nu)
    for(no=nu:Nopt)
        i=no-nu+1; 
        H(no,nu)=h(i);
    end
end
% obliczamy macierz WR
WR=inv(H'*H+rDu*eye(Nu))*H'; 