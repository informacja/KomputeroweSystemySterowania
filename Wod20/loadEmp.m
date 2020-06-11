function [Xtime, YempTemperatura] = loadEmp(Ystart, Ystop, filename)

%     Ystart = 2017; % rok od którego zacz¹æ wczytwyaæ temperatury
%     Ystop = 2019; % do którego roku w³¹cznie
%     filename = 'Raport_Pomiarow_anonim.xls';

    if ~exist('x', 'var')% jeœli jeszcze nie s¹ w workspace, za³aduj dane z pliku
        x = []; y = [];

        for i = Ystart:Ystop
            %     x = importfile('Raport_Pomiarow_anonim.xls', string(i));
            x = [x; importRaport_Pomiarow(filename, string(i))];
            y = [y; size(x(:, 7), 1)];
        end

    end

    %% table to matrix

    id_lp   = find(string(x.Properties.VariableNames) == 'Lp');
    id_date = find(string(x.Properties.VariableNames) == 'Datapoboru');
    id_celc = find(string(x.Properties.VariableNames) == 'Wynik');
    id_place = find(string(x.Properties.VariableNames) == 'Miejscepoboru');

    c       = string(x{:, id_celc}); % Temperatura w Celcjuszach
    place   = strtrim(strip(string(x{:, id_place}), 'both')); % Miejsce
    lp      = x{:, id_lp}; % Liczba porzadkowa
    d       = x{:, id_date}; % Data poboru probki
    date    = datetime(d, 'InputFormat', 'dd.MM.yyyy  ', 'Format', 'yyyy-MM-dd');

    exp = '(\d*[.,]*\d*)';
    celc = regexp(c, exp, 'match');

    if length(celc) == length(date)
        disp('Test - D³ugoœæ danych wejœciowyjch jest zgodna [ok] ');
    else
        length(celc)
        length(data)
        error('Ro¿ne dlugosci tablic, dane mog¹ byæ niepoprawne [Error]')
     
    end

    %%

    z = [];
    Tdate = [];

    for i = 1:size(celc, 1)
        if contains(place(i), 'SUW £ukanowice')
            if double(celc{i}) > 30 | double(celc{i}) < 1 % czy próbka spe³nia relane  kryteria   fprintf( "Próbka %s nie spe³nia kryterium, temperatura wody %.1f *C\ns",  date(i), double(celc{i}))
                  fprintf( 'Próbka %s nie spe³nia kryterium, temperatura wody %.1f *C\ns',  date(i), double(celc{i}))
            else
                Tdate = [Tdate; date(i)];
                d = double(celc{i});
                z = [z; d];
            end
        end
        %     z = x(i,1)
        %     c(i)
    end
%     z(size(Tdate, 1)) = 0;
length(Tdate);
Tdate = Tdate(1:size(z,1));
length(Tdate);
% figure(3),plot( Tdate, z ); title('Surowe dane');

YempTemperatura = z;
Xtime = Tdate;
end
% size (z)
% size (Tdate)

% plot( Tdate, z ); title("Surowe dane");
% figure, plot( Tdate, detrend(z,0) ); title("detrend mean");
% figure, plot( Tdate, detrend(z,1) ); title("detrend linear");
% figure, plot( Tdate, detrend(z,2) ); title("detrend quadratic trend");
