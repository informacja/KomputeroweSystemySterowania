% Cel:
% Czy cykl roczny jest istotny? 7 14 21
%   CSV wygenerowaæ
%   naprawiæ wykresy
%   dobraæ harmoniczne
clear all; %Cwicz2Regr

if (verLessThan('matlab', '8.4')) % for MATLAB 2010
%     filename = '2018_space.csv'; delimiterIn = ',';
%     data = importdata(filename, delimiterIn)'; % wektor pionowy (szafa)
    load('time_data.mat')
    x = datenum(time'); tR = [1:max(x) - min(x) + 2]; % time Real
    L = length(x);
    T = max(x) - min(x) + 1;
    L = length(data); x = [1:L]; tR = x; % generowanie numeru próbki dla MatLab 2010
    time = x;

else
    [time, data] = loadEmp(2013, 2019, 'Raport_Pomiarow_anonim.xls'); % lata od 2013 do 2019
    x = datenum(time'); tR = [1:max(x) - min(x) + 2]; % time Real

    if (length(x) == length(unique(x))) fprintf(1, "ok")% To do
    else length(unique(x)); error ("kilka pomiarów jednego dnia");
    end
end
%----------------------------------------------------------------------


%---------
% srednia
%---------

% ------------------------------------------------------------------------
Ldanych = length(x); % size(x, 2);
Yemp = dtrend(data, 1); sigYf = std(Yemp);
figure(1), subplot(2, 1, 1), plot(time, data); axis('tight'); 
ylabel('Temperatura *C'); legend(sprintf('Œrednia = %i', mean(Yemp))); title('Dziedzina czasu');
Ah = abs(fft(Yemp / Ldanych)); Ah = Ah(1:round(Ldanych / 2)); %nie dzia³a round w matlab 2010
nrOm = find(Ah > 2.5 * mean(Ah)) - 1;
subplot(2, 1, 2), plot(Ah); axis('tight'); title('Dziedzina czêstotliwoœci');

% G³êbokie myœlenie = prezentacja danych + wyobra¿enia
kol1 = 'r*'; if(Ldanych > 80) kol1 = 'r.'; end
figure(2), subplot(2, 1, 1);
plot(x, Yemp, kol1); xlabel('Numer wêz³a'); ylabel('odchy³ temperatury [*C]'); title('Trend usuniêty'); hold on

% z = input(' ? jaka to funkcja ? !!!  <Ent> - co mogloby byc ?');

% za³. funkcji harmonicznej
% Nrom = nrOm; %[1 2 5 20 36 42 134 236 500 600];

om = 2 * pi / T * nrOm;
% %% Projekt modelu - oblicz FId
[FId, Ldanych, Lh, Kd] = modelRharm(x, om); % za³. funkcji harmonicznej
% =============== wybieramy model: Lhm i Km ===================
Lhm = 3;
% =============================================================
Km = 2 * Lhm; % Przyjety (arbitralnie) rzad modelu  Kolumny Macierzy
FI(:, 1:Km) = FId(:, 1:Km);
% dalej ju¿ tylko numeryka i grafika
Gd = FI' * FI;
G = inv(Gd);
Ao = G * FI' * Yemp; % Wsp.A=inv(FItransp*FI)*FItransp*Yempiryczne
% sprawdz modelu
Yo = FI * Ao;
Em = Yemp - Yo;

if (1) E = Em; LdE = Ldanych;
else E(x') = Em; LdE = length(E); % uzupe³niamy braki zerami;
end

varE = (E' * E) / (LdE - Km);
sigEo = sqrt(varE);

% ============= model AR reszt ======================
alf = E(2:end)' * E(1:end - 1) / ((LdE - 1) .* varE);
Talf = -1 / log(alf); txalf = '\alpha';
v(LdE, 1) = 0;
v(1, 1) = 0; for(n = 2:LdE) v(n, 1) = E(n - 1) * alf; end;
sigZ = std(E - v);
% ........................................................
% obliczanie sigYo i sigYe;
KA = G * varE; % macierz kowariancji wspolcz.
Lduzych = 0;

for (n = 1:Ldanych)
    varYo = FI(n, :) * KA * FI(n, :)';
    sigYo(n, 1) = sqrt(varYo);
    sigYe(n, 1) = sqrt(varYo + varE);
    if (abs(E(n)) > sigYe(n)) Lduzych = Lduzych + 1; end
end

Kz = sigEo * sqrt(1 - alf^2);

uLd = Lduzych / Ldanych * 100;
% ................................
hold on; if(Ldanych > 30) kol = 'b'; else kol = 'bo-'; end
plot(x, Yemp, 'k.', x, Yo, kol, x, E, 'r', x, E - v, 'g', x, v, 'b', [x(1) x(end)], [0 0], 'k:'); %axis('tight');
%xlabel(sprintf('Ldanych=%d sigYf=%.3f sigEo=%.3f', Ldanych, sigYf, sigEo));
plot(x, Yo + sigYo, 'b--', x, Yo - sigYo, 'b--');
plot(x, Yo + sigYe, 'm--', x, Yo - sigYe, 'm--');
hold off; axis('tight')
xlabel(sprintf('Ldanych=%d sigYf=%.3f sigEo=%.3f udz.DyzychE=%.1f%% %s=%.3f T_%s=%.3f', Ldanych, sigYf, sigEo, uLd, txalf, alf, txalf, Talf));
% ......................................................
fprintf(1, '\nWspolczynniki i ich statystyki tSt');
m = 0;

for (k = 1:Kd)

    if (k > Km) tSt(k) = 0; A(k) = 0;
    else m = m + 1; tSt(k) = abs(Ao(m)) / KA(m, m); A(k) = Ao(m);
    end

    % fprintf(1, '\nA(%-2d)=%-9.3f tSt=%-7.3f ....... Af=%-7.3f', k, A(k), tSt(k), Af(k)); % DO ZATWIERDZENIA
    fprintf(1, '\nA(%-2d)=%-9.3f tSt=%-7.3f', k, A(k), tSt(k));
end

fprintf(1, '\n alf=%.4f Talf=%.3f Kz=%.3f\n', alf, Talf, Kz)

xmin = min(x); xmax = max(x);
Ldim = 2000; Xmin = xmin; Xmax = xmin + T * 1.3;
dv = (Xmax - Xmin) / (Ldim - 1); v = [Xmin:dv:Xmax];
z0 = 0; z(Ldim) = 0;

for (i = 1:Ldim)

    fi = modelRharm(v(i), om(1:Lhm)); % za³. funkcji harmonicznej

    yo(i, 1) = fi * Ao;
    varYv = fi * KA * fi';
    sigYv(i, 1) = sqrt(varYv);
    sigYE(i, 1) = sqrt(varYv + varE);
    z(i) = alf * z0 + Kz * randn; % z - E symulowane
    z0 = z(i);
end

figure(2), subplot(2, 1, 2), plot(x, Yemp, kol1, v, yo, 'r', x, E, 'k', v, z, 'g'); hold on;
%input(' co dalej ?? ');
plot(v, yo + sigYv, 'b:', v, yo - sigYv, 'b:');
plot(v, yo + sigYE, 'm:', v, yo - sigYE, 'm:');
xlabel(sprintf('Ldanych=%d sigYf=%.3f sigEo=%.3f Km=%d udz.DyzychE=%.1f%%', Ldanych, sigYf, sigEo, Km, uLd));
axis('tight'); title('AR Poza polem korelacji'); hold off;

% end
return
