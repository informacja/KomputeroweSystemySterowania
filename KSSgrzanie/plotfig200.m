%plotfig200;
nsh=shift*nDt;
tSk=ceil((nskok-d1+1)/nDt)*DtR; %tSko=tSk-d+1;
figure(200); subplot(3,1,1);
xdy=[1: length(yd)]*DtR; xdu=xdy;
plot(xdy,yd,'k.-',xdu,Fid(:,2),'g.:',xdy,ytf,'m.'); axis('tight');
ax=axis; hold on;
plot([tSk tSk], ax(3:4),'b:',[tSk tSk]+nsh, ax(3:4),'b-.'); hold off;
axis([0.99*ax(1) 1.01*ax(2) ax(3:4)]);

subplot(3,1,2);
%xdy=[d+(przyr+1+shift)*nDt:nDt:dyf]*dt; xdu=[d1+(przyr+shift)*nDt:nDt:nuf-nDt]*dt;
%yd=Y(d:nDt:dyf)'-y0;
%yp=Y(d-nDt:nDt:dyf-nDt)'-y0; ud=U(d1:nDt:nuf-nDt)'-u0;
dp=(przyr)*nDt;
xdy=[d+dp:nDt:d+dp+(LFi-1)*nDt]*dt;  %xdy=[d+dp:nDt:dyf]*dt;
xdu=[d1+dp:nDt:nuf-nDt]*dt;
ny0=[d1+dp:nDt:d+dp+nDt];

plot(ny0(1:end-1)*dt,Y(ny0(1:end-1)),'g.-',xdy,yd,'k.-',xyt,ytm,'c.-',xyt,yt,'b.:',...
    xdy(1:end)-dt,yn,'k.-',xdu(kk:end),Fid(kk:end,2),'r:',xdy,ytf,'m.');
axis('tight'); ax=axis;
hold on;
plot([tSkok tSkok], ax(3:4),'b--',[d d]*dt,ax(3:4),'g:',[d+nDt+nsh d+nDt+nsh]*dt,ax(3:4),'b-.',...
    [dyf dyf]*dt,ax(3:4),'b:',[d1 d1]*dt,ax(3:4),'r:',[nuf-nDt nuf-nDt]*dt,ax(3:4),'r:',ax(1:2),[0 0],'k:');
axis([0.99*ax(1) 1.01*ax(2) ax(3:4)]); hold off;

%delr=nd*dt;
alf=A(1,1); Tor=-DtR/log(alf);
if(jestK==0) Kor=A(2,1)/(1-alf); end;
xlabel(sprintf('ndf=%d J=%.3g Ryy=%.3g RCN=%.3g delr=%.3f Kor=%.2f Tor=%.2f dr=%d',ndf,J(ndf),Kyy(ndf),RCN(ndf),delr,Kor,Tor,d));
lw=3; lk=3; lfp=(lw-1)*lk;
if(kf==1) lk=2; lfp=(lw-1)*lk;
else subplot(lw,lk,lfp+3); plot(RCNx,'k'); axis('tight'); xlabel('RCNx');
end
%if(kk==1) subplot(3,1,3); plot(xdy,E,'k'); axis('tight'); xlabel('yd-ytf');
%else
subplot(lw,lk,lfp+1); plot(xdy(2:end),E,'k',xyt,Em,'r'); xlabel('(yd-ytf) k; (y-yt) r');
axis('tight'); ax=axis; hold on; plot(ax(1:2),[0 0],'k:'); hold off;
if(Koniec==0)
    xtJ=([d01:nDt:d]-2*nDt)*dt;
    subplot(lw,lk,lfp+2); plot(xtJ,JE(1:ndf),'k',xtJ,Jm(1:ndf),'r',round(delmin/nDt)*DtR,J(ndmin),'mo'); axis('tight'); 
    xlabel(sprintf('Jx(delR): delmin=%d',delmin));
end
%subplot(lk,lw,lfp+4); plot(Kyux,'k'); axis('tight'); xlabel('Kyux');
%end
axx=1;
%if(abs(delr-0.4)<1.e-4) Ao=A; fio=fi; yo=y;  end
%input(' ?? ');
PokazId=PokazId-1; 
if(PokazId==0)
    PokazId=input(' Co ile kroków pokaz ?? Wpisz <Ent> - do nastêpnego; w=0 - do konca; w>0 - po w.krokach; w<0 - DEBUG i po |w|.krokach ?? ');
    if(isempty(PokazId)) PokazId=1; end
    if(PokazId<0) PokazId=abs(PokazId); dbstop('nic');  nic; end; %nic(coile);
end


