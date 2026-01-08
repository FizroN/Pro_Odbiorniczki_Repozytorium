clear all; close all; clc;

% =========================================================================
% 1. KONFIGURACJA
% =========================================================================
APP_MODE = 1;              % 1 = Nagrywanie, 2 = Odczyt pliku .wav
BANDWIDTH_MHz = 1.024;     % Dostępne: 1.024, 5, 20

if BANDWIDTH_MHz == 20
    Fs = 20e6; RECORD_DURATION = 0.1; Gain = 40;
elseif BANDWIDTH_MHz == 5
    Fs = 5e6; RECORD_DURATION = 5.0; Gain = 40;
elseif BANDWIDTH_MHz == 1.024
    Fs = 1024000; RECORD_DURATION = 50; Gain = 73;
end

Fc = 112.6e6;  
DataFolder = ''; % Bieżący folder

% Tworzenie nazwy pliku bezpośrednio z rozszerzeniem .wav
FileName = sprintf('Airplane_IQ_%.1fMHz_%.1fs_%.0fHz_Krakow.wav', BANDWIDTH_MHz, RECORD_DURATION, Fc);
FilePath = fullfile(DataFolder, FileName);

% =========================================================================
% 2. POBIERANIE / WCZYTYWANIE DANYCH
% =========================================================================

if APP_MODE == 1
    % --- NAGRYWANIE ---
    % Zakładamy, że funkcja zwraca sygnał do zmiennej receivedSignal
    % Jeśli funkcja wymaga ścieżki do zapisu .mat, przesyłamy tymczasową lub modyfikujemy ją
    [receivedSignal, metadata] = RecordPlutoSignal('', Fc, Fs, RECORD_DURATION, Gain);
    
    % Przygotowanie danych do zapisu (I + Q jako Stereo)
    if size(receivedSignal, 2) > size(receivedSignal, 1), receivedSignal = receivedSignal.'; end
    wavData = [real(receivedSignal), imag(receivedSignal)];
    
    % Normalizacja do zakresu [-1, 1] dla formatu .wav
    maxAmp = max(abs(wavData(:)));
    if maxAmp > 0, wavData = wavData / maxAmp; end
    
    % Zapis bezpośrednio do .wav
    try
        audiowrite(FilePath, wavData, Fs);
        fprintf('--- NAGRYWANIE I ZAPIS .WAV ZAKOŃCZONY ---\n');
        fprintf('Plik: %s\n', FilePath);
    catch ME
        error('Błąd zapisu pliku .wav: %s', ME.message);
    end

elseif APP_MODE == 2
    % --- ODCZYT ---
    if isfile(FilePath)
        SelectedFile = FilePath;
    else
        [fName, fPath] = uigetfile('*.wav', 'Wybierz plik sygnału IQ (.wav)');
        if isequal(fName,0), error('Anulowano.'); end
        SelectedFile = fullfile(fPath, fName);
    end
end

% =========================================================================
% 3. WERYFIKACJA I DEMODULACJA
% =========================================================================
% Nie używaj clear all tutaj!

% 1. Rekonstrukcja IQ z poprawionym znakiem
[raw, fs] = audioread(FilePath);
x = raw(:,1) + 1i*raw(:,2); % PLUS jest kluczowy

% 2. Automatyczne szukanie częstotliwości środkowej (Piku stacji)
Nfft = 2^18; 
[Pxx, f_axis] = periodogram(x, kaiser(length(x),10), Nfft, fs, 'centered');
[~, idx] = max(Pxx);
fc_detected = f_axis(idx);

fprintf('Wykryto najsilniejszy sygnał na offsecie: %.2f kHz\n', fc_detected/1e3);

% Ustawienie fc do demodulacji na to, co faktycznie znaleźliśmy
fc = fc_detected; 
demod = 3;

% Wyświetlenie widma, żebyś widział czy masz piki
figure(10);
plot(f_axis/1e3, 10*log10(Pxx));
hold on; plot(fc/1e3, 10*log10(Pxx(idx)), 'ro'); % Zaznaczenie piku
grid on; xlabel('Częstotliwość [kHz]'); ylabel('Moc [dB]');
title('Widmo nagranego sygnału IQ');

