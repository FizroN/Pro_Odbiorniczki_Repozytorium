% Program_Nadajnik.m
% Wysyła ramkę danych (Preambuła + Tekst) w pętli przez ADALM-PLUTO

clear all; close all;
run Ustawienia.m; % Wczytaj wspólne parametry

%% 1. Przygotowanie Danych (Ramka)
% Definicja konstelacji
[IQcodes, Nstates, Nbits, ~] = IQdef(modtype); %

% Kodowanie tekstu na numery symboli
numbers_data = text2numbers(tekst_testowy, Nbits); %

% Generowanie preambuły
[preamble_nums, ~] = modtype2header(modtype); %

% Łączenie preambuły i danych
numbers_tx = [preamble_nums; numbers_data];
N_syms = length(numbers_tx);

% Konwersja numerów na symbole IQ
IQk_tx = numbers2IQ(numbers_tx, modtype, IQcodes); %

% Wstawianie zer (Upsampling)
IQ0 = zeros(1, N_syms * K);
IQ0(1:K:end) = IQk_tx;

% Formowanie impulsów (Filtr TX)
IQn_tx = conv(IQ0, K * hpsf); %

% Normalizacja sygnału, aby uniknąć przesterowania
IQn_tx = IQn_tx / max(abs(IQn_tx)) * 0.5;
txWaveform = transpose(IQn_tx); % Nadajnik Pluto wymaga wektora kolumnowego

fprintf('Nadajnik uruchomiony. Konfiguracja ADALM-PLUTO...\n');

%% 2. Konfiguracja Nadajnika Pluto
tx = sdrtx('Pluto'); % Podłącz do Pluto
tx.RadioID = 'usb:0'; % Użyj domyślnego radia
tx.CenterFrequency = fcar;
tx.BasebandSampleRate = fs;
tx.Gain = 0; % Wzmocnienie w dB (0 jest zwykle bezpieczne)

% Wysyłaj sygnał w pętli
transmitRepeat(tx, txWaveform);

fprintf('Transmisja rozpoczęta na %d MHz. Naciśnij Ctrl+C, aby zatrzymać.\n', fcar/1e6);
fprintf('Wysyłanie tekstu: "%s"\n', tekst_testowy);

% Pętla, aby skrypt się nie zakończył (zatrzymaj przez Ctrl+C)
while true
    pause(1);
end