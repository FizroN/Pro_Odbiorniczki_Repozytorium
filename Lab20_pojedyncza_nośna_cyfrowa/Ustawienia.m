% Ustawienia.m
% Wspólne parametry dla nadajnika, odbiornika i loopbacku

%% Parametry Modulacji
modtype = 'QPSK'; % Wybierz: '2PAM', '4PAM', 'QPSK', '16QAM' itd.
tekst_testowy = 'Hello World! 0123456789';

%% Parametry Filtrowania i Próbkowania
r = 0.35;                   % Współczynnik "roll-off" filtra
K = 8;                      % Liczba próbek na symbol (SPS)
Ns = 10;                    % Długość filtra (w symbolach)
Npsf = Ns * K + 1;          % Długość filtra w próbkach
Mpsf = (Npsf - 1) / 2;      % Połowa długości filtra
rctype = 'sqrt';            % Typ filtra: 'sqrt' dla TX i RX

%% Parametry Systemu
fs = 400e3;                 % Szybkość próbkowania (Hz) - 400 kHz
fcut = fs / (2 * K);        % Częstotliwość odcięcia filtra
fcar = 100e6;               % Częstotliwość nośnej (Hz) - 100 MHz (dla Pluto)

%% Generowanie Filtra
% Projektujemy filtr SRRC (pierwiastkowy)
hpsf = firrcos(Npsf - 1, fcut, r, fs, 'rolloff', rctype); %