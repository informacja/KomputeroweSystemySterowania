% model.m
function Ym=model(K,T,d,U,n,Yp)
global dt;
if(n<=d) u=U(1); else u=U(n-d); end
dY=-Yp/T+K/T*u; 
Ym=Yp+dY*dt; 
