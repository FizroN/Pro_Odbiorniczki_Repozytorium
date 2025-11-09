% Program_Odbiornik.m
% Odbiera i dekoduje dane na żywo z ADALM-PLUTO

clear all; close all;
run Ustawienia.m; % Wczytaj wspólne parametry

fprintf('Odbiornik uruchomiony. Konfiguracja ADALM-PLUTO...\n');

%% 1. Przygotowanie Odbiornika
% Definicja konstelacji
[IQcodes, Nstates, Nbits, ~] = IQdef(modtype); %

% Wygeneruj preambułę, której będziemy szukać
[preamble_nums, ~] = modtype2header(modtype); %
preamble_IQk = numbers2IQ(preamble_nums, modtype, IQcodes); %
len_preamble = length(preamble_IQk);

% Przygotuj filtr dopasowany (jest już w 'hpsf' z Ustawienia.m)
matchedFilter = hpsf;

% Przygotuj obiekty do wizualizacji (opcjonalnie, ale pomocne)
cd = comm.ConstellationDiagram('Title', 'Odebrane Symbole', ...
    'ReferenceConstellation', IQcodes);

%% 2. Konfiguracja Odbiornika Pluto
rx = sdrrx('Pluto');
rx.RadioID = 'usb:0';
rx.CenterFrequency = fcar;
rx.BasebandSampleRate = fs;
rx.SamplesPerFrame = 50000; % Rozmiar bufora
rx.OutputDataType = 'double';
rx.GainSource = 'AGC Slow Attack'; % Automatyczna kontrola wzmocnienia

% Zmienne stanu pętli
bufor_nadmiarowy = complex(zeros(Npsf-1, 1));
znaleziono_ramke = false;
licznik_symboli = 0;
dataBuffer = [];

fprintf('Rozpoczynanie nasłuchu na %d MHz...\n', fcar/1e6);

%% 3. Główna Pętla Odbiorcza
while true
    
    % Krok 1: Odbiór bufora danych
    [bufor_rx, ~] = rx();
    if isempty(bufor_rx)
        continue;
    end
    
    % Krok 2: Filtr Dopasowany
    % Łączymy z nadmiarowymi próbkami z poprzedniej pętli, aby
    % nie zgubić symboli na granicy buforów
    rx_samples = [bufor_nadmiarowy; bufor_rx];
    IQn_filtered = conv(rx_samples, matchedFilter);
    bufor_nadmiarowy = rx_samples(end-Npsf+2:end); % Zapisz nadmiar na następną pętlę
    
    % Krok 3: Synchronizacja (Wykrywanie Preambuły)
    % To jest uproszczona wersja. Szukamy korelacji z preambułą
    % po *downsamplingu* (co K-tą próbkę).
    
    % Downsampling (prosta wersja - bierzemy co K-tą próbkę)
    % Lepsza metoda wymagałaby pętli synchronizacji symbolu (np. Gardnera)
    % Zaczynamy od opóźnienia filtrów
    IQk_downsampled = IQn_filtered(Mpsf*2+1 : K : end);
    
    % Korelacja krzyżowa, aby znaleźć preambułę
    [korelacja, lagi] = xcorr(IQk_downsampled, preamble_IQk);
    [max_corr, idx] = max(abs(korelacja));
    
    prog_detekcji = length(preamble_IQk) * 0.5; % Próg detekcji
    
    if max_corr > prog_detekcji
        % ZNALEZIONO SYGNAŁ!
        
        % Znajdź początek ramki
        start_idx = lagi(idx) + 1;
        
        % Użyj preambuły do oszacowania błędu fazy/CFO
        % (Uproszczenie: tylko błąd fazy)
        odebrana_preambula = IQk_downsampled(start_idx : start_idx + len_preamble - 1);
        
        % Estymacja błędu fazy
        phase_error = angle(sum(odebrana_preambula .* conj(preamble_IQk.')));
        
        % Korekcja całego strumienia (uproszczenie)
        IQk_corrected = IQk_downsampled * exp(-1j * phase_error);
        
        % Wizualizuj skorygowane symbole
        step(cd, IQk_corrected);
        
        % Ekstrakcja danych (symbole *po* preambule)
        data_start_idx = start_idx + len_preamble;
        dane_symbole = IQk_corrected(data_start_idx : end);
        
        % Dekodowanie
        numbers_rx = IQ2numbers(dane_symbole, modtype); %
        
        % Konwersja na tekst
        % Musimy wiedzieć, ile symboli tekstu wysłano (z Ustawienia.m)
        len_data = length(text2numbers(tekst_testowy, Nbits)); %
        
        if length(numbers_rx) >= len_data
            numbers_data_rx = numbers_rx(1:len_data);
            text_out = numbers2text(numbers_data_rx, Nbits); %
            
            % Wyświetl wynik
            fprintf('----------------------------------------\n');
            fprintf('ZNALEZIONO RAMKĘ! Odebrano: %s\n', text_out);
            fprintf('----------------------------------------\n');
        end
    end
end