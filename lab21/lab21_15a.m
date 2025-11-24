% lab21_ex_FMradio_encoder_stereo_RDS.m
  clear all; close all;
  
  ifigs = 1;       % 0/1 figures No/Yes
  FMradio_params   % read parameters from the included file

% Read two monophonic audio signals - audioread or wavread
if(1)
  [x1, fx1 ] = audioread('GOODBYE.WAV');  x1=x1.'; % plot(x1);
  [x2, fx2 ] = audioread('DANKE.WAV');    x2=x2.'; % soundsc(x1,fx1);
else % for test only
  Nsec = 5; fx1 = faudio; fx2 = fx1; N = Nsec*fx1; t = 1/fx1*(0:N-1);
  x1 = 0.5*cos( 2*pi*(5000*t - 4900/(2*pi*0.2)*cos(2*pi*0.2*t)));
  x2 = 0.5*cos( 2*pi*(5000*t + 4900/(2*pi*0.2)*cos(2*pi*0.2*t)));
% x2 = zeros(1,length(t)); % for testing cross-talk between channels
end

% Re-sampling to FM radio audio frequency - faudio
  [N1,D1] = rat(faudio/fx1,1e-6); x1 = resample(x1,N1,D1);   % resampling
  [N2,D2] = rat(faudio/fx2,1e-6); x2 = resample(x2,N2,D2);   %
  Nx1=length(x1); Nx2=length(x2); Nx = max( Nx1, Nx2 );
  x1 = [ x1 zeros(1,Nx-Nx1) ];                               % appending zeros
  x2 = [ x2 zeros(1,Nx-Nx2) ];                               % 
% x2 = zeros(1,Nx); % for testing cross-talk between channels in the receiver

% Filters were designed for faudio=25000 Hz, our signals can have different frequency
% Pre-emphasis, flat freq response to 2.1 kHz, than increasing 20 dB per decade
% x1 = filter(b_pre,a_pre,x1); x2 = filter(b_pre,a_pre,x2); 

% Up-sampling to FM radio service frequency - fs
  [N,D] = rat(fs/faudio,1e-6); x1 = resample(x1,N,D); x2 = resample(x2,N,D);
   
% RDS BITS #########################################
% --- ZADANIE 21.15 CZĘŚĆ 1: Nadawanie własnego tekstu ---
% A. Generation of RDS bits from text
  
  % 1. Twój tekst
  my_text = 'HALO HALO TU RADIO GEMINI '; 
  
  % 2. Konwersja tekstu na ciąg bitów (ASCII -> Binary)
  % dec2bin tworzy macierz znaków '0' i '1'
  bin_matrix = dec2bin(my_text, 8); 
  
  % 3. Przekształcenie w jeden długi wektor liczb 0 i 1
  % Odejmujemy 48, bo w ASCII znak '0' ma wartość 48, a '1' ma 49.
  % Transpozycja (') jest ważna, żeby czytać wierszami (znak po znaku).
  bits_vector = double(bin_matrix') - 48;
  bits_vector = bits_vector(:)'; % Upewniamy się, że to wektor poziomy
  
  % 4. Powielanie tekstu, żeby wypełnić cały czas trwania sygnału
  Nx = length(x1); 
  Nsps = fs/fsymb; 
  Nrds_needed = ceil(Nx/Nsps); % Tyle bitów potrzebujemy łącznie
  
  % Powielamy wektor bitów tyle razy, ile trzeba
  num_repeats = ceil(Nrds_needed / length(bits_vector));
  rds = repmat(bits_vector, 1, num_repeats);
  
  % Przycinamy do dokładnej długości
  rds = rds(1:Nrds_needed);

% --------------------------------------------------------
% B. Differential encoding of RDS bits
  Nrds = length(rds); r = 0;
  for k=2:Nrds, r(k)=abs(rds(k)-r(k-1));
  end
% C. Generation of bi-phase impulses: -1 -> +1, +1 -> -1, 
  bip = [-1 1; 1 -1 ]; brds = bip( r+1,: ); brds=brds'; brds=brds(:);
% D. Zero insertion (appending) between -1/+1 impulses
  Nr = ceil(fs/fsymb *length(r)); uprds = zeros(1,Nr); t = dt*(0:Nr-1);
  clk = abs(diff(sign(sin(2*pi/Tsymb*t))))/2; clk = [clk 0];
  indx=1;
  for n=1:Nr
      if(clk(n)==1) uprds(n)=brds(indx); indx=indx+1;
      end
  end
% E. Signal smoothing using pulse shaping filter (PSF), synchro with pilot   
  prds = conv( uprds, hpsf); prds = prds(Npsf:Npsf+Nx-1);
  clear t clk r brds uprds;
% RDS BITS #########################################

% Generation of multiplex (MPX) signal of FM radio
  n=0:Nx-1; alpha = 2*pi*fpilot/fs*n;
  x = 0.9*(x1+x2)+0.1*cos(alpha)+0.9*(x1-x2).*cos(2*alpha)+0.05*prds.*cos(3*alpha);
  A = max(x); x=x/A;

  if(ifigs==1)
     m=1:20*210;
     figure; plot(m,cos(alpha(m)/16),'r-',m,prds(m),'b-'); grid; title('prds(n)'); pause
     m=1:2000; N = length(x); f = fs/N*(0:N-1);
     figure, plot(x(m)); title('MPX(n) signal'); xlabel('n'); grid; pause;
     figure; plot(f,20*log10(abs(fft(x)))); xlabel('f (Hz)'); title('FFT of MPX signal'); grid; pause
     figure, spectrogram(x, 512, 512-64, 512, fs); title('STFT of MPX signal'); pause
  end

% Final frequency modulation of the complex carrier in the base-band
  BW = 160000;                % overall FM radio bandwidth > 2*fmax=2*60000 Hz
  fmax = 60000;               % maximum modulating frequency
  df = (BW/(2*fmax)-1)*fmax,  % from Carson's bandwidth rule
  beta = df/fmax,             % modulation index
  fc = 0;                     % carrier frequency of FM radio in the base-band
  x = exp( j*2*pi*( fc/fs*n + df*cumsum(x)/fs ) ); x=x.';  % FM modulation
  I = real(x); Q = imag(x);                                % I, Q components
  % --- KONIEC ENKODERA (MODYFIKACJA: Wyraźny zapis) ---

% TU JEST KLUCZOWY MOMENT ZAPISU
% Zapisujemy składowe I oraz Q jako dwa kanały pliku WAV
output_filename = 'Test_Radio_Hello.wav';
fprintf('Zapisywanie sygnału do pliku: %s ... ', output_filename);

audiowrite(output_filename, [I Q], fs); % [cite: 2062]

fprintf('GOTOWE!\n');
fprintf('Teraz uruchom dekoder, aby odczytać ten plik.\n');
  if(ifigs==1)
     m=1:5000; N = length(x); f = fs/N*(0:N-1);
     figure, plot(m,I(m),m,Q(m)); title('FM Signal Re/Im = IQ'); xlabel('n'); grid; pause;
     figure; plot(f,20*log10(abs(fft(x)))); xlabel('f (Hz)'); title('FFT of FM Signal'); grid; pause
     figure, spectrogram(x, 512, 512-64, 512, fs); title('STFT of FM Signal'); pause
  end
