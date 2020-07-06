%  Cel:
% Czy cykl roczny jest istotny? 7 14 21
% +  CSV wygenerowaæ - > time_data.mat
% +  naprawiæ wykresy
%   dobraæ harmoniczne (cel czy cykl roczny jest istotny) brak
%   harmonicznych
%   uzup³niæ brak pomiarów z regresji do fft
%    montecarlo
clear all; %Cwicz2Regr



% if (verLessThan('matlab', '8.4')) % for MATLAB 2010
% %     filename = '2018_space.csv'; delimiterIn = ',';
% %     data = importdata(filename, delimiterIn)'; % wektor pionowy (szafa)
%     load('time_data.mat')
%     x = datenum(time'); tR = [1:max(x) - min(x) + 2]; % time Real
%     L = length(x);
%     T = max(x) - min(x) + 1;
%     L = length(data); x = [1:L]; tR = x; % generowanie numeru próbki dla MatLab 2010
% %     time = x;
% 
% else % for MATLAB 2020
%     [time, data] = loadEmp(2013, 2019, 'Raport_Pomiarow_anonim.xls'); % lata od 2013 do 2019
%     x = datenum(time'); tR = [1:max(x) - min(x) + 2]; % time Real
% 
%     if (length(x) == length(unique(x))) fprintf(1, 'ok - 1 pomiar dziennie')% To do
%     else length(unique(x)); error ('kilka pomiarów jednego dnia lub inny powód powtórzeñ w wektorze [time]');  end
% end
%----------------------------------------------------------------------

%[rTemp,Temp,nTim,czas,lTim,lrT]=dane1()
[Yemp,data,time,x,lTim,lrT]=dane(); 
% [size(Yemp);size(data);size(time);size(x);size(lTim);size(lrT)] % debug
% ------------------------------------------------------------------------
% UWAGA do widma furerowiskiego musi byæ spe³nione za³o¿enie o równomiernym
% próbkowaniu, poni¿sza linia jest tylko gdzy harmoniczne ju¿ znamy 
% Yemp = data; x = time; % podmian na dane z dziurami do modelu regresjnego 
Ldanych = length(x); % size(x, 2);
T = max(x) - min(x) + 1; sred_temp = mean(Yemp);
Yemp = dtrend(Yemp, 1); sredAfterDentrend = mean(Yemp);
sigYf = std(Yemp);
figure(1), subplot(2, 1, 1), plot(time, data); axis('tight'); 
ylabel('Temperatura *C'); xlabel(sprintf("Œrednia temperatura %d", sred_temp)); title('Dziedzina czasu');
Ah = abs(fft(Yemp / Ldanych));  Ah = Ah(1:round(Ldanych / 6)); %nie dzia³a round w matlab 2010
acept_level = 4.5 * mean(Ah); A(1:length(Ah)) = acept_level;
f0 = find(max(Ah) == Ah); % cykl roczny
harm = [Ah(f0:f0:end)] % To DO uwzglêdniaj harmoniczne roczne w modelu
nrOm = find(Ah > acept_level) - 1; 
subplot(2, 1, 2),   plot([1:length(Ah)],Ah(1:end),'k.-',nrOm+1,Ah(nrOm+1),'g^',1:length(Ah),A,'r--');  

% plot(Ah([7 14]),'go'); plot([acept_level acept_level],'g:');
% axis('tight');
ylabel('Amplituda'); xlabel(sprintf(' Proponowane %i silne harmoniczne jako bazowa nr = %s ',length(nrOm),  mat2str(nrOm))); title('Dziedzina czêstotliwoœci'); 
legend('g^ Proponowane\nro Wybrane');
%---------
% srednia
%---------
% input('srednia ruchoma');
%% sr +-2
% [m,n] = size(Yemp);
% y1 = Yemp';
% y2 = zeros(m,n);
% v = 2;
% for i = v+1:(n-v)    
%     sum = y1(i);
%     for (w = 1:v)
%        sum = sum + (y1(i-w) + y1(i+w));
%     end
%     y2(i) = sum/(v+3);
% end

% Yemp = y2


% G³êbokie myœlenie = prezentacja danych + wyobra¿enia
kol1 = 'r*'; if(Ldanych > 80) kol1 = 'r.'; end
figure(2), subplot(2, 1, 1);
plot(x, Yemp, kol1);  ylabel('[*C]'); title('G³êbokie myœlenie = prezentacja danych + wyobra¿enia'); hold on

% z = input(' ? jaka to funkcja ? !!!  <Ent> - co mogloby byc ?');

% za³. funkcji harmonicznej
nrOm = [7 14 21 28 56]; %nrOm; %[1 2 5 20 36 42 134 236 500 600];

om = 2 * pi / T * nrOm;
% %% Projekt modelu - oblicz FId
[FId, Ldanych, Lh, Kd] = modelRharm(x, om); % za³. funkcji harmonicznej
% dodaæ dane z dziurami
% =============== wybieramy model: Lhm i Km ===================
Lhm = 5;
% =============================================================
Yemp=Yemp';
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
axis('tight'); title('co dalej ??'); hold off;

% end
return
