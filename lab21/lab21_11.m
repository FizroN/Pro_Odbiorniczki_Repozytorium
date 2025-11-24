% --- ZADANIE 21.11: Wizualizacja filtrów radia FM ---

% 1. Wczytanie parametrów i projektowanie filtrów
% Uruchamiamy skrypt dostarczony z książką (Listing 21.8)
% Jeśli nie masz pliku 'FMradio_params.m', daj znać - trzeba go będzie stworzyć.
try
    FMradio_params; 
catch
    error('Nie znaleziono skryptu FMradio_params.m! Upewnij się, że jest w folderze.');
end

% Lista filtrów do zwizualizowania (nazwy zmiennych z FMradio_params)
% hpsf      - Filtr kształtujący (SRRC) dla RDS
% hLPaudio  - Filtr dolnoprzepustowy audio
% hBP19     - Filtr pilota 19 kHz
% hBP38     - Filtr podnośnej stereo 38 kHz
% hBP57     - Filtr RDS 57 kHz
% b_de, a_de - Filtr de-emfazy

figure('Name', 'Analiza filtrów Radia FM');

% --- 1. Filtr De-emfazy (IIR) ---
subplot(3,2,1);
[H, W] = freqz(b_de, a_de, 1024, fs);
plot(W, 20*log10(abs(H)));
title('De-emfaza (IIR)'); xlabel('Hz'); ylabel('dB'); grid on;
xlim([0 15000]); % Skupiamy się na paśmie audio

% --- 2. Filtr Audio LP (FIR) ---
subplot(3,2,2);
[H, W] = freqz(hLPaudio, 1, 1024, fs);
plot(W, 20*log10(abs(H)));
title('Audio Low-Pass'); xlabel('Hz'); grid on;
xlim([0 20000]);

% --- 3. Filtr Pilota 19kHz ---
subplot(3,2,3);
[H, W] = freqz(hBP19, 1, 4096, fs);
plot(W, 20*log10(abs(H)));
title('Pilot 19 kHz BP'); xlabel('Hz'); grid on;
xlim([15000 23000]); % Zoom na okolice 19k

% --- 4. Filtr Stereo 38kHz ---
subplot(3,2,4);
[H, W] = freqz(hBP38, 1, 4096, fs);
plot(W, 20*log10(abs(H)));
title('Stereo 38 kHz BP'); xlabel('Hz'); grid on;
xlim([34000 42000]); % Zoom na okolice 38k

% --- 5. Filtr RDS 57kHz ---
subplot(3,2,5);
[H, W] = freqz(hBP57, 1, 4096, fs);
plot(W, 20*log10(abs(H)));
title('RDS 57 kHz BP'); xlabel('Hz'); grid on;
xlim([53000 61000]); % Zoom na okolice 57k

% --- 6. Filtr Kształtujący SRRC (dla RDS) ---
subplot(3,2,6);
[H, W] = freqz(hpsf, 1, 1024, fs);
plot(W, 20*log10(abs(H)));
title('RDS Pulse Shaping (SRRC)'); xlabel('Hz'); grid on;