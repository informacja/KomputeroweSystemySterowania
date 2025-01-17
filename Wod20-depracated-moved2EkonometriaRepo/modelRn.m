function [E,varZ,sigZo,Yo,G,FI,Kd,Ao] = modelR(x,Yemp,Km,Ldanych,Lh,om,k) %Nrom sigZf,Lhm,Kd,T,
%% Projekt modelu - oblicz FId
% za�. funkcji harmonicznej
% Nrom = [1 2 5 20 36 42 134 236 500 600];

for (nh = 1:Lh)
    k = k + 1; FId(:, k) = sin(om(nh) * x'); k = k + 1; FId(:, k) = cos(om(nh) * x');
end

Kd = k;
% wybieramy model
FI(:, 1:Km) = FId(:, 1:Km);
% dalej ju� tylko numeryka i grafika
Gd = FI' * FI;
G = inv(Gd);
Ao = G * FI' * Yemp; % Wsp.A=inv(FItransp*FI)*FItransp*Yempiryczne
% sprawdz modelu
Yo = FI * Ao;
E = Yemp - Yo; varZ = (E' * E) / (Ldanych - Km);
sigZo = sqrt(varZ);
end