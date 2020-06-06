% plik RegDMC.m
function [DU]=RegDMC(Yref,WR,h,yn,Nopt,dU, Ur, DelR, Umax, Umin, DUreg,n)
global DtR jestF111 typReg; %nPrdYref
lh=length(h); ldU=length(dU); lU=length(Ur); Up=Ur(lU-1); 
% Liczymy y swob. 
n0=n-lh; if(n0<=0) n0=1; end
ynp=h(lh)*Ur(n0); % yn biez¹ca wartosc; ynp - bedzie teoretyczna wart. bie¿¹ca
for(i=2:lh-1)
    ndU=ldU-i-DelR+2;
    if(ndU<1) break; end
    ynp=ynp+h(i)*dU(ndU);
end
vpr=yn-ynp; % prognoza ZOH
% Liczymy y swob. 
n=ldU; np=DelR; nk=Nopt+DelR-1; 
np1=np-1;
%n0=n+np-lh; if(n0<=0) n0=1; end
for(p=np:nk)
    ypl(p,1)=vpr; 
% Uwaga !!! tu h(1)=0; h(2:lh) odp.skok. bez opóŸnienia
% ypl(p)=Sum(h(i=2:lh)*dU(n+p-[p+1:p+1+lh]))
    ih1=p+2-np;
    %for(i=p+1:lh) 
    %    iu=n+p-i;
    for(i=2:lh-1) 
        iu=n+p-i-np1; 
        if(iu>=n) continue; end; % Tu dU=dUpl=0;
        if(iu<1) iu=1; break; end
        ypl(p,1)=ypl(p,1)+h(i)*dU(iu); 
    end
    ypl(p,1)=ypl(p,1)+h(lh)*Ur(iu); 
end
ldYr=length(Yref(:,1)); Yref(ldYr:nk,1)=Yref(ldYr,1); 
DU=-WR*(ypl(np:nk,1)-Yref(np:nk,1));
if(Up+DU(1,1)>Umax) DU(1,1)=Umax-Up; end
if(Up+DU(1,1)<Umin) DU(1,1)=Up-Umin; end
if(abs(DU(1,1))>DUreg) DU(1,1)=sign(DU(1,1))*DUreg; end
if(mod(n,30)==0 && jestF111) 
    figure(111); %jestF111=1;
    for(i=1:length(DU(:,1))) Ur(n-1+i)=Ur(n-2+i)+DU(i,1); end
    Ur(n-1+i:n-1+Nopt+DelR)=Ur(n-1+i); 
    %tR=[n:n+Nopt+DelR]*DtR;
    %plot(tR,Ur(n-1:n-2+Nopt+DelR),'g',[n:n+Nopt+DelR-1]*DtR,Yref(1:Nopt+DelR-1,1),'m',[np:nk]*DtR,ypl(np:nk,1),'k'); 
    plot([n-1:n-2+Nopt+DelR]*DtR,Ur(n-1:n-2+Nopt+DelR),'g',...
        [n:n+Nopt+DelR-2]*DtR,Yref(1:Nopt+DelR-1,1),'m',...
        (n-1+[np:nk])*DtR,ypl(np:nk,1),'k',...
        (n-1+[np:nk])*DtR,ypl(np:nk,1)-vpr,'k:',...
        [n n]*DtR,[yn yn],'bo',[n n]*DtR,[ynp ynp],'r*');    
    axis('tight'); title(typReg);
    hold on; ax=axis; 
    plot([n n]*DtR,ax(3:4),'k:',[n+DelR n+DelR]*DtR,ax(3:4),'r:')
    aa=1; 
end