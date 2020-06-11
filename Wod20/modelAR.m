function [E,varZ,sigZo,Yo,G,FI,Kd,Kz,alf] = model(x,Nrom,Yemp,sigZf,Lhm,Km,Ldanych,Lh,Kd,T,om,k);
%% Projekt modelu - oblicz FId
% za³. funkcji harmonicznej
% Nrom = [1 2 5 20 36 42 134 236 500 600];

for (nh = 1:Lh)
    k = k + 1; FId(:, k) = sin(om(nh) * x'); k = k + 1; FId(:, k) = cos(om(nh) * x');
end

Kd = k;
% wybieramy model
FI(:, 1:Km) = FId(:, 1:Km);
% dalej ju¿ tylko numeryka i grafika
Gd = FI' * FI;
G = inv(Gd);
Ao = G * FI' * Yemp; % Wsp.A=inv(FItransp*FI)*FItransp*Yempiryczne
% sprawdz modelu
Yo = FI * Ao;
E = Yemp - Yo; varZ = (E' * E) / (Ldanych - Km);
sigZo = sqrt(varZ);



% ============= model AR reszt ======================
alf = E(2:end) * E(1:end - 1)' / ((Ldanych - 1) .* varZ);
% alf = E(2:end) * E(1:end - 1); DO ZATWIERDZENIA transpozycja
v(Ldanych) = 0;
v(1) = 0; for(n = 2:Ldanych)
E(n - 1) .* alf';
end;
Kz = std(E - v);
hold on; if(Ldanych > 30) kol = 'k.'; else kol = 'ko'; end
plot(x, Yo, kol); %axis('tight');
xlabel(sprintf('Ldanych=%d sigZf=%.3f sigZo=%.3f', Ldanych, sigZf, sigZo));
hold off;
% obliczanie sigYo i sigYe;
KA = G * varZ; % macierz kowariancji wspolcz.
Lduzych = 0;

for (n = 1:Ldanych)
    varYo = FI(n, :) * KA * FI(n, :)';
    sigYo(n, 1) = sqrt(varYo);
    sigYe(n, 1) = sqrt(varYo + varZ);
    if (abs(E(n)) > sigYe(n)) Lduzych = Lduzych + 1; end
end

uLd = Lduzych / Ldanych * 100;
hold on;
plot(x, Yo + sigYo, 'b--', x, Yo - sigYo, 'b--');
plot(x, Yo + sigYe, 'm--', x, Yo - sigYe, 'm--');
hold off; axis('tight')
xlabel(sprintf('Ldanych=%d sigZf=%.3f sigZo=%.3f udz.DyzychE=%.1f%%', Ldanych, sigZf, sigZo, uLd));
fprintf(1, '\nWspolczynniki i ich statystyki tSt');

for (k = 1:Kd)

    if (k > Km)tSt(k) = 0; A(k) = 0;
    else tSt(k) = abs(Ao(k)) / KA(k, k); A(k) = Ao(k);
    end

    % fprintf(1, '\nA(%-2d)=%-9.3f tSt=%-7.3f ....... Af=%-7.3f', k, A(k), tSt(k), Af(k)); % DO ZATWIERDZENIA
    fprintf(1, '\nA(%-2d)=%-9.3f tSt=%-7.3f ....... Af=%-7.3f', k, A(k), tSt(k), alf(k));
end





xmin = min(x); xmax = max(x);
Ldim = 2000; Xmin = xmin; Xmax = xmax * 1.3;
dv = (Xmax - Xmin) / (Ldim - 1); v = [Xmin:dv:Xmax]; z0 = 0;
sgZ = Kz; Kz = sigZo * sqrt(1 - alf^2);

for (i = 1:Ldim)
    k = 0;

    for (nh = 1:Lhm)
        k = k + 1; fi(k) = sin(om(nh) * v(i));
        k = k + 1; fi(k) = cos(om(nh) * v(i));
    end

    yo(i, 1) = fi * Ao;
    varYv = fi * KA * fi';
    sigYv(i, 1) = sqrt(varYv);
    sigYE(i, 1) = sqrt(varYv + varZ);
    z(i) = alf * z0 + Kz * randn; z0 = z(i); % DO ZATWIERDZENIA niem mo¿na przypisaæ macierzy do indexu
end

plot(x, Yemp, kol1, v, yo, 'r', v, z, 'g.'); hold on;
%input(' co dalej ?? ');
plot(v, yo + sigYv, 'b:', v, yo - sigYv, 'b:');
plot(v, yo + sigYE, 'm:', v, yo - sigYE, 'm:');
xlabel(sprintf('Ldanych=%d sigZf=%.3f sigZo=%.3f Km=%d udz.DyzychE=%.1f%%', Ldanych, sigZf, sigZo, Km, uLd));
axis('tight'); hold off;

end