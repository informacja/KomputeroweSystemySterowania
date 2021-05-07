function [u,du,ep,Sp_epp]=RegPID(y,Yref,Kr,TI,TD,Dt,Umin,Umax,DU,ep,Sp_epp,U0,Up,jaki)
e=Yref-y(3); 
if(jaki=='p')
    Sp_epp=Sp_epp+e; 
    u=U0+Kr*(e+Dt/TI*Sp_epp+TD/Dt*(e-ep)); 
    ep=e; 
else
    u=U0+Kr*(e-ep+Dt/TI*e+TD/Dt*(e-2*ep+Sp_epp)); 
    Sp_epp=ep; ep=e;
end
if(u<Umin) u=Umin; else if(u>Umax) u=Umax; end, end
du=(u-Up); if(abs(du)>DU) du=sign(du)*DU; u=Up+du; end