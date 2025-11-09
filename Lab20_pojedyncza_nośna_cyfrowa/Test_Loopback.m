% Test_Loopback.m
% Weryfikacja całego łańcucha DSP (TX -> Kanał -> RX) w symulacji

clear all; close all;
run Ustawienia.m; % Wczytaj wspólne parametry

fprintf('Uruchomiono test pętli zwrotnej (Loopback)...\n');
fprintf('Wiadomość do wysłania: %s\n', tekst_testowy);

%% 1. NADAJNIK (TX)
% Definicja konstelacji
[IQcodes, Nstates, Nbits, ~] = IQdef(modtype); %

% Kodowanie tekstu na numery symboli
numbers_data = text2numbers(tekst_testowy, Nbits); %

% Generowanie preambuły
[preamble_nums, ~] = modtype2header(modtype); %
fprintf('Długość preambuły: %d symboli\n', length(preamble_nums));
fprintf('Długość danych: %d symboli\n', length(numbers_data));

% Łączenie preambuły i danych
numbers_tx = [preamble_nums; numbers_data];
N_syms = length(numbers_tx);

% Konwersja numerów na symbole IQ
IQk_tx = numbers2IQ(numbers_tx, modtype, IQcodes); %

% Wstawianie zer (Upsampling)
IQ0 = zeros(1, N_syms * K); %
IQ0(1:K:end) = IQk_tx; %

% Formowanie impulsów (Filtr TX)
IQn_tx = conv(IQ0, K * hpsf); %

%% 2. KANAŁ (Symulacja)
% Używamy funkcji 'IQdisturb' do symulacji kanału
% Ustawmy realistyczne, ale trudne parametry
SNR = 20;           % Stosunek sygnał-szum w dB
chanG = 1.0;        % Wzmocnienie kanału
chanPh = pi/8;      % Przesunięcie fazy kanału (45 stopni)
carDF = 0.001;      % Błąd częstotliwości nośnej (CFO) jako ułamek fs
carDPh = pi/3;      % Błąd fazy nośnej
ADCdt = 0.2;        % Błąd taktowania ADC (ułamek okresu próbkowania)
Nskip = Mpsf;       % Pominięcie próbek na brzegach

% Zastosuj zakłócenia
IQn_channel = IQdisturb(IQn_tx, SNR, chanG, chanPh, carDF, carDPh, ADCdt, Nskip); %

%% 3. ODBIORNIK (RX)
% Filtr dopasowany (Matched Filter)
IQn_rx = conv(IQn_channel, hpsf); %

% --- Synchronizacja (w symulacji) ---
% W prawdziwym odbiorniku to najtrudniejsza część.
% Tutaj "oszukujemy", bo znamy opóźnienia z filtrów.
% Całkowite opóźnienie = filtr TX (Mpsf) + filtr RX (Mpsf)
total_delay = 2 * Mpsf;

% Ekstrakcja symboli (Downsampling)
% Musimy zacząć od pierwszej próbki po opóźnieniu filtrów
IQk_rx = IQn_rx(total_delay + 1 : K : end); %

% --- Korekcja (w symulacji) ---
% W prawdziwym RX potrzebowalibyśmy tu pętli PLL do korekcji
% fazy i częstotliwości. W symulacji z 'IQdisturb', błędy
% carDF i carDPh nie są idealnie skompensowane, ale spróbujmy.

% Przycinamy do oryginalnej długości (jeśli nadmiarowe próbki)
IQk_rx = IQk_rx(1:N_syms);

% Dekodowanie symboli IQ na numery
% UWAGA: 'IQ2numbers' nie koryguje błędów fazy/CFO!
% Działa tylko dla idealnych danych lub PAM. Dla QPSK/16QAM
% ten test może się nie powieść bez pętli korekcyjnej.
numbers_rx = IQ2numbers(IQk_rx, modtype); %

% Konwersja numerów na tekst (pomijamy preambułę)
% Znajdź długość preambuły i odetnij ją
len_preamble = length(preamble_nums);
numbers_data_rx = numbers_rx(len_preamble + 1 : end);

text_out = numbers2text(numbers_data_rx, Nbits); %

%% 4. Weryfikacja
fprintf('----------------------------------------\n');
fprintf('Odebrana wiadomość: %s\n', text_out);

if strcmp(strtrim(tekst_testowy), strtrim(text_out))
    fprintf('SUKCES: Wiadomości są identyczne!\n');
else
    fprintf('BŁĄD: Wiadomości są różne.\n');
end

% Wizualizacja konstelacji (po korekcji opóźnienia, przed decyzją)
figure;
plot(real(IQk_rx), imag(IQk_rx), 'b.');
hold on;
plot(real(IQcodes), imag(IQcodes), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
title('Odebrana konstelacja (przed decyzją)');
xlabel('I'); ylabel('Q'); grid on;
legend('Odebrane symbole', 'Idealna konstelacja');