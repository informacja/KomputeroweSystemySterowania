function [Yemp,Yteor,sigZf,Af]=obiekt(x,sigZf) %,alfa,v0)
%Af=[1,-0.8,1.2,0.3,0,0,0,0 0 0 0 0 0 0]; ld=length(x);
ld=length(x(:,1)); Lx=length(x(1,:));
%Af=[1, 150, 0.0, -15, -8.0, 0.0, 0.5e1, 0 0 0 0 0 0 0];
Af=[1,-0.8,1.2,0.3,0,0,0,0 0 0 0 0 0 0]; ld=length(x);
%Af=[1,-6.0,3.2,0.3,0,0,0,0 0 0 0 0 0 0];
kf=length(Af); Kf=kf;
fi(1:ld,1)=1;
for(k=2:Kf) fi(:,k)=x.^(k-1); end
% Liczymy wyjscie
Yteor=zeros(ld,1);
for k=1:Kf
    Yteor=Yteor+Af(k)*fi(:,k);
end
if sigZf<0
    sigZf=-sigZf;
    Dy=max(Yteor)-min(Yteor);
    sigZf=Dy*sigZf; % odchylenie stand Z;
end
% Obliczamy zaklocenia
v=sigZf*randn(ld,1);
Yemp=Yteor+v;
return;
% =================================




if(nargin<3 || alfa<1.e-4) 
    Hammer=0; 
    if(Lx==1)
        %Af=[1, 150, 0.0, -15, -8.0, 0.0, 0.5e1, 0 0 0 0 0 0 0]; 
        Af=[1,-0.8,1.2,0.3,0,0,0,0 0 0 0 0 0 0]; ld=length(x);
        %Af=[1,-6.0,3.2,0.3,0,0,0,0 0 0 0 0 0 0]; 
        kf=length(Af); Kf=kf; 
        fi(1:ld,1)=1; 
        for(k=2:Kf) fi(:,k)=x.^(k-1); end
    else
        Af=zeros(1,15); 
        fi(1:ld,1)=1; k=1;
        for(i=1:2) k=k+1; fi(:,k)=x(:,i); end % czlony 1.go rzedu
        Af(1)=1; Af(2)=-0.8; Af(3)=22; % k=2 i 3
        k=k+1; fi(:,k)=x(:,1).*x(:,2); % k=4; 
        Af(4)=10; 
        for(i=1:2) k=k+1; fi(:,k)=(x(:,i).^2); end % czl.2go rzedu
        Af(5)=0; Af(6)=0; 
%        k=k+1; fi(:,k)=(x(:,1).^2).*x(:,2); % mieszane 3.go rzedu
%        Af(7)=0; Af(8)=0;    
%        k=k+1; fi(:,k)=x(:,1).*x(:,2).^2; % miesz.3go rzedu
%        Af(9)=0; Af(10)=0;    
%        for(i=1:2) k=k+1; fi(:,k)=(x(:,i).^3); end % czl.3go rzedu
        Kf=k; 
    end   
else
    Hammer=1; 
    Af=[1, 15, alfa]; kf=length(Af); Kf=kf-1; 
end
Yteor=zeros(ld,1);
for k=1:Kf
   Yteor=Yteor+Af(k)*fi(:,k); 
end
if sigZf<0
    sigZf=-sigZf;
    Dy=max(Yteor)-min(Yteor);
    sigZf=Dy*sigZf;
end
% Obliczamy zaklocenia
if(Hammer) 
    K=sigZf*sqrt(1-alfa^2);
    if(nargin==3 || v0==0.) 
        v0=0; 
        for(n=1:100) v0=alfa*v0+K*randn; end, 
    end, % rozbiegowy
    for(n=1:ld) v0=alfa*v0+K*randn; v(n,1)=v0; end, 
else v=sigZf*randn(ld,1); 
end
Yemp=Yteor+v;