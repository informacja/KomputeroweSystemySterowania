function [mi_sr miw etaw]=lepkw(Twg,Twz,Tp)
global ro0 wExp aa bb Bs As juzByleta 
Nargin=3; %Nargin=nargin; %
Twz=40; Twg=70;
if(juzByleta==0)
    % Dane wg http://fizyk.ifpk.pk.edu.pl/tabele/wLepWody.htm
    %juzByleta=1;
    Aetaw=[20 10.02 	35 	7.20 	50 	5.48 	65 	4.36 ...
        21 	9.78 	36 	7.05 	51 	5.39 	66 	4.30 ...
        22 	9.55 	37 	6.94 	52 	5.31 	67 	4.23 ...
        23 	9.33 	38 	6.80 	53 	5.25 	68 	4.18 ...
        24 	9.09 	39 	6.68 	54 	5.14 	69 	4.12 ...
        25 	8.91 	40 	6.54 	55 	5.06 	70 	4.06 ...
        26 	8.68 	41 	6.42 	56 	4.98 	71 	4.00 ...
        27 	8.51 	42 	6.30 	57 	4.90 	72 	3.95 ...
        28 	8.22 	43 	6.18 	58 	4.82 	73 	3.89 ...
        29 	8.15 	44 	6.06 	59 	4.75 	74 	3.85 ...
        30 	7.98 	45 	5.97 	60 	4.68 	75 	3.80 ...
        31 	7.81 	46 	5.86 	61 	4.60 ...
        32 	7.63 	47 	5.76 	62 	4.54 ...
        33 	7.48 	48 	5.66 	63 	4.47 ...
        34 	7.33 	49 	5.56 	64 	4.42];
    Tw=Aetaw(1:2:end); Etae=Aetaw(2:2:end);
    
    %Tz=Twg+(Twg-Tp)*exp(-Kkal/(Fw*cw));
    Mie=Etae/ro0.*(1+wExp*(Tw-20));
    aa=polyfit(Tw,Etae,3);
    bb=polyfit(Tw,Mie,3);
    As=[aa(1)/4 aa(2)/3 aa(3)/2 aa(4)];
    Bs=[bb(1)/4 bb(2)/3 bb(3)/2 bb(4)];
end
clear etaw miw mi1s;
if(juzByleta==0)
    juzByleta=1;
    eta=aa(4)+aa(3)*Tw+aa(2)*Tw.^2+aa(1)*Tw.^3; 
    mi=bb(4)+bb(3)*Tw+bb(2)*Tw.^2+bb(1)*Tw.^3; 
    figure(200); subplot(1,2,1); plot(Tw,eta,'k.',Tw,Etae,'ro');
    subplot(1,2,2); plot(Tw,mi,'k.',Tw,Mie,'ro');
end
if(nargout>1)
    etaw(1)=aa(4)+aa(3)*Twz+aa(2)*Twz.^2+aa(1)*Twz.^3;
    miw(1)=bb(4)+bb(3)*Twz+bb(2)*Twz.^2+bb(1)*Twz.^3;
end
mi1s=Twz*(Bs(4)+Bs(3)*Twz+Bs(2)*Twz^2+Bs(1)*Twz^3);
mi_sr=mi1s; 
if(Nargin>1)
    if(nargout>1)
        etaw(2)=aa(4)+aa(3)*Twg+aa(2)*Twg.^2+aa(1)*Twg.^3;
        miw(2)=bb(4)+bb(3)*Twg+bb(2)*Twg.^2+bb(1)*Twg.^3;
    end
    mi2s=Twg*(Bs(4)+Bs(3)*Twg+Bs(2)*Twg.^2+Bs(1)*Twg.^3);
    mi_sr=(mi1s-mi2s)/(Twz-Twg);
end

