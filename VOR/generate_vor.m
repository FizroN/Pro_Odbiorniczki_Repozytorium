clear all; close all; clc;

% =========================================================================
% KONFIGURACJA GENERATORA VOR
% =========================================================================
Fs = 100000;          % Częstotliwość próbkowania (100 kHz - jak w analizatorze)
Duration = 10;        % Czas trwania (sekundy)
Azimuth_deg = 90;     % <--- TUTAJ USTAW KĄT (np. 90 stopni = Wschód)
NoiseLevel = 0.05;    % Poziom szumu (0.0 = idealny, 0.1 = realny)

% =========================================================================
% MATEMATYKA SYGNAŁU
% =========================================================================
t = 0:1/Fs:Duration-1/Fs; % Oś czasu

% 1. Sygnał Zmienny (Variable) - 30 Hz AM
% W VOR ten sygnał jest przesunięty w fazie o kąt azymutu
Azimuth_rad = deg2rad(Azimuth_deg);
variable_signal = 0.3 * cos(2*pi*30*t - Azimuth_rad); 

% 2. Sygnał Odniesienia (Reference) - 30 Hz FM na podnośnej 9960 Hz
% FM: nośna 9960 Hz, sygnał modulujący 30 Hz (kosinus), indeks modulacji 16
reference_signal = 0.3 * cos(2*pi*9960*t + 16 * sin(2*pi*30*t));

% 3. Identyfikacja Morse'a (1020 Hz) - "TEST" (- . ... -)
% Uproszczona wersja (ciągły ton przerywany kluczowaniem)
morse_tone = 0.1 * cos(2*pi*1020*t);
% Kluczowanie (prosta maska: włącz/wyłącz co jakiś czas)
keying = double(square(2*pi*1*t) > 0); % Bip co 0.5 sekundy
ident_signal = morse_tone .* keying;

% 4. Składanie całości (Sygnał Baseband - wokół 0 Hz)
% Dodajemy stałą składową (DC) dla AM
signal_clean = (1 + variable_signal) .* (1 + reference_signal + ident_signal);
sig_complex = signal_clean + 1j*hilbert(signal_clean)*0.0; 

% !!! NOWOŚĆ: Przesuwamy sygnał o 25 kHz w prawo (symulacja Offset Tuningu) !!!
OffsetFreq = 25000; % 25 kHz
sig_complex = sig_complex .* exp(1j*2*pi*OffsetFreq*t); 

% (Dalej tak samo - dodanie szumu)
noise = NoiseLevel * (randn(size(t)) + 1j*randn(size(t)));
receivedSignal = sig_complex + noise;

% Normalizacja
receivedSignal = receivedSignal / max(abs(receivedSignal));

% =========================================================================
% ZAPIS DO PLIKU
% =========================================================================
FileName = sprintf('VOR_Synthetic_Az%d.wav', Azimuth_deg);
OutputFolder = 'Signals';
if ~exist(OutputFolder, 'dir'), mkdir(OutputFolder); end
FullPath = fullfile(OutputFolder, FileName);

% Zapis jako WAV (Stereo: L=Real, R=Imag)
wavData = [real(receivedSignal(:)), imag(receivedSignal(:))];
audiowrite(FullPath, wavData, Fs);

fprintf('--- GENERATOR ZAKOŃCZONY ---\n');
fprintf('Wygenerowano sygnał dla azymutu: %d stopni\n', Azimuth_deg);
fprintf('Zapisano plik: %s\n', FullPath);
fprintf('Teraz użyj tego pliku w swoim analizatorze!\n');