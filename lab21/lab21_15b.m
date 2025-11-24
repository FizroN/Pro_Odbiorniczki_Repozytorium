% lab21_ex_FMradio_decoder_stereo_RDS.m
  clear all; close all;
 
  FMradio_params   % reading parameters from the included file

% READ IQ SIGNAL - synthesized or recorded
  % --- POCZĄTEK DEKODERA ---
% READ IQ SIGNAL - Wczytujemy plik wygenerowany przez Enkoder

input_filename = 'Test_Radio_Hello.wav';

if exist(input_filename, 'file')
    fprintf('Wczytywanie pliku: %s ...\n', input_filename);
    [y, fs] = audioread(input_filename); % [cite: 2185]
else
    error('Nie znaleziono pliku %s! Najpierw uruchom Enkoder.', input_filename);
end                    % synth
% [y,fs] = audioread('SDRSharp_FMRadio_96000kHz_IQ_one.wav',[1,1*250000]); % record

% FM DEMODULATION
y = y(:,1) - sqrt(-1)*y(:,2);            % IQ --> complex
dy = (y(2:end).*conj(y(1:end-1)));       % calculation of instantaneous frequency
y = atan2(imag(dy), real(dy)); clear dy; % frequency demodulation

% DECODING MONO L+R SIGNAL
ym = filter( hLPaudio, 1, y );           % low-pass filtration of L+R (mono) signal
ym = ym(1:fs/faudio:end);                % leaving only every fs/faudio-th sample
disp('LISTENING TO MONO: Left + Right'); soundsc(ym,faudio); pause
w = ym-mean(ym); w=w/max(abs(w)); audiowrite('FM_mono.wav',w,faudio); clear w;

% CARRIER RECOVERY
  p = filter(hBP19,1,y);                            % extracting 19 kHz pilot
  theta = zeros(1,length(p)+1);                     %# 
  omega = theta; omega(1) = 2*pi*fpilot/fs;         %# double PLL on 19 kHz
  mi1 = 0.0025; mi2 = mi1^2/4; p=p/max(abs(p));     %# see chapter on AM
  for n = 1 : length(p)                             %#
      pherr = -p(n)*sin(theta(n));                  %# 
      theta(n+1) = theta(n) + omega(n) + mi1*pherr; %# 
      omega(n+1) = omega(n) + mi2*pherr;            %#
  end                                               %#
  c1(:,1)    = cos(theta(1:end-1)/16);              %@ for RDS synchro
  c1PSF(:,1) = cos(theta(1:end-1)/16+phasePSF);     %@ 19 kHz/16 = 1187.5 Hz 
  c38(:,1) = cos(2*theta(1:end-1));                 % L-R carrier 38 kHz
  c57(:,1) = cos(3*theta(1:end-1));                 % RDS carrier 57 kHz
  clear p; clear theta; clear freq;
        
% DECODING STEREO
  if(1) ys = filter(hBP38,1,y); delay = (L/2)/(fs/faudio); % extraction of L-R signal
  else ys=y; delay=0; end                     %(optional BP filtration around 38 kHz)
  ys = real(ys .* c38); clear c38;            % L-R signal: 38kHz --> 0kHz 
  ys = filter( hLPaudio, 1, ys );             % low-pass filtration
  ys = ys(1:fs/faudio:end);                   % leaving every fs/faudio-th sample
  ym = ym(1:end-delay); ys=2*ys(1+delay:end); % synchronization of L+R and L-R
  clear ymm yss; 
  y1 = 0.5*( ym + ys ); y2 = 0.5*( ym - ys ); clear ym ys;    % recovering L and R
  % y1 = filter(b_de,a_de,y1); y2 = filter(b_de,a_de,y2);     % de-emphasis
  disp('LISTENING TO STEREO'); soundsc([y1'; y2'],faudio); pause
  maxi = max( max(abs(y1)),max(abs(y2)) );
  audiowrite('FM_Radio_stereo.wav',[ y1/maxi y2/maxi ],faudio); clear y1 y2;

% DECODING RDS
% A. Initial operations
  if(1) y = filter(hBP57,1,y); c1 = c1PSF; end  % extraction of RDS (optional BP filter)
  y = y .* c57; clear c57;                      % frequency conversion: 57kHz-->0kHz
  y = filter( hpsf, 1, y );                     % low-pass filer - SRRC pulse shaping 
        n = 20000 : 20000 + 50*210-1;
        figure; plot(n,c1(n),'r-',n,y(n),'b-'); grid; title('c1(n), y(n)'); pause
% B. Signal correlation with clock c1 shifted 16x by 1/fpilot (exact timing recovery) 
  nstart = 20000; Nmany = 50*210;               % skipping time of PLL loop adaptation
  for n = 1 : round(fs/fpilot) : round(fs/fsymb)+1;
      synchro(n)=sum((y(nstart+n-1:nstart+n-1+Nmany-1).* c1(nstart:nstart+Nmany-1)));  
  end
  [v, indx] = max(synchro);                     % finding maximum  
  y = y(indx:end); c1=c1(1:end-indx+1);         % final signal synchronization
        n = 20000 : 20000 + 50*210-1;
        figure; plot(n,c1(n),'r-',n,y(n),'b-'); grid; title('Matched c1(n), y(n)'); pause
% C. Signal multiplication with synchronized carrier - coherent detection
  y = y .* c1; y = y/max(y);  
% D. Moving average of last fs/Tsymb samples
  h = ones(1,round(fs/fsymb))/round(fs/fsymb); y = filter(h,1,y); 
% E. Finding negative and positive slopes of synchronized clock
  N = length(c1); M=5; M2=(M-1)/2; maxi1 = []; maxi2 = [];
  for n = M2+1:N-M2
    if( (c1(n-2)>0) && (c1(n-1)>0) && (c1(n+1)<0) && (c1(n+2)<0) ) % negative slope
        maxi1 = [ maxi1 n ];
    end    
    if( (c1(n-2)<0) && (c1(n-1)<0) && (c1(n+1)>0) && (c1(n+2)>0) ) % positive slope
       maxi2 = [ maxi2 n ];
    end
  end
  if( std(y(maxi1)) > std(y(maxi2)) ) maxi=maxi1; else maxi=maxi2; end % best fitting?
        figure; n=1:length(y); 
        plot(maxi1,c1(maxi1),'bo',maxi2,c1(maxi2),'ro',n,c1,'r-',n,y,'b-'); title('Detect'); pause
% F. Changing signal levels: {-1, +1} --> {0,1} 
  rdsbits = (-sign( y( maxi ) )+1)/2; clear c1 maxi maxi1 maxi2;
% G. Differential decoding with taking into account bi-phase signal nature
  rdsbits = abs(rdsbits(2:end) - rdsbits(1:end-1)); rdsbits = rdsbits(2:2:end)
% H. Storing detected bits to disk for further analysis by outside program
  % --- ZADANIE 21.15 CZĘŚĆ 2: Odczyt wiadomości RDS ---
% Podmieniamy SCRX na rdsbits, bo tak nazywa się zmienna w Twoim pliku.

fprintf('\n--- ODCZYTANA WIADOMOŚĆ RDS ---\n');

% 1. Upewniamy się, że bity są binarne (0/1)
bits = rdsbits;
bits(bits < 0.5) = 0;
bits(bits >= 0.5) = 1;

% 2. Szukamy tekstu, sprawdzając różne przesunięcia bitowe
% (nie wiemy, gdzie zaczyna się bajt, więc sprawdzamy 8 możliwości)

for shift = 0:7
    % Bierzemy bity z przesunięciem
    shifted_bits = bits(1+shift : end);
    
    % Obcinamy do pełnych bajtów
    n_bytes = floor(length(shifted_bits) / 8);
    if n_bytes == 0, continue; end
    
    bits_cut = shifted_bits(1 : n_bytes*8);
    
    % Formujemy macierz: 8 wierszy x n_bytes kolumn
    bytes_mat = reshape(bits_cut, 8, n_bytes);
    
    % Konwersja bin -> dec (ASCII)
    % Wagi bitów: najstarszy bit (128) na górze macierzy
    powers = [128 64 32 16 8 4 2 1];
    ascii_codes = powers * bytes_mat;
    
    % Zamiana na znaki
    decoded_str = char(ascii_codes);
    
    % Wyświetlamy fragment dla każdego przesunięcia
    % Żeby nie zaśmiecać, wyświetlamy tylko jeśli tekst wygląda na czytelny
    % (np. zawiera litery i spacje, a nie kody sterujące < 32)
    readable_chars = sum(ascii_codes >= 32 & ascii_codes <= 126);
    score = readable_chars / length(ascii_codes);
    
    fprintf('Offset %d: "%s" (Czytelność: %.1f%%)\n', shift, decoded_str(1:min(60,end)), score*100);
end

fprintf('\n---------------------------------\n');   
