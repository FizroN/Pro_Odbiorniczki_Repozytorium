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
% A. Generation of dummy RDS bits for encoding
  Nx = length(x1); Nsps=fs/fsymb; Nrds = ceil(Nx/Nsps); rds = [];
  for k=1:ceil(Nrds/10), rds = [ rds 1 1 1 1 0 0 0 1 1 0 ];
  end
% rds = round(rand(1,Nrds)); % for test
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
  clear t clk;
% RDS BITS #########################################

% Generation of multiplex (MPX) signal of FM radio
  n=0:Nx-1; alpha = 2*pi*fpilot/fs*n;

% --- ZADANIE 21.12: Wizualizacja etapów generacji RDS (POPRAWIONY) ---
% Wklej ten kod w pliku 'lab21_ex_FMradio_encoder_stereo_RDS.m'
% DOKŁADNIE PRZED linią tworzącą zmienną 'x' (zazwyczaj x = 0.9 * ...)

figure('Name', 'Generacja sygnału RDS (Analiza)');

% 1. Surowe bity po kodowaniu różnicowym: r(k)
subplot(2,2,1);
% Używamy length(r) zamiast samego end wewnątrz min()
L_r = length(r);
stem(r(1:min(20, L_r)), 'filled'); 
title('r(k) - Bity różnicowe'); 
xlabel('k (indeks bitu)'); grid on;

% 2. Sekwencja Bi-phase (Manchester): brds(k)
subplot(2,2,2);
temp_brds = brds(:); % Spłaszczamy wektor
L_brds = length(temp_brds);
stem(temp_brds(1:min(40, L_brds)), 'filled');
title('brds(k) - Kodowanie Bi-phase'); 
xlabel('k (próbki symboli)'); grid on;

% 3. Sygnał po wstawieniu zer (Upsampling): uprds(n)
subplot(2,2,3);
L_uprds = length(uprds);
stem(uprds(1:min(2000, L_uprds)), 'filled', 'MarkerSize', 2);
title('uprds(n) - Wstawione zera'); 
xlabel('n (próbki sygnału)'); grid on;
ylim([-1.2 1.2]); 

% 4. Sygnał po filtracji kształtującej: prds(n)
subplot(2,2,4);
L_prds = length(prds);
plot(prds(1:min(2000, L_prds)), 'LineWidth', 1.5);
title('prds(n) - Po filtrze SRRC'); 
xlabel('n (próbki sygnału)'); grid on;
xlim([0 2000]);

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
  audiowrite('FMRadio_IQ_250kHz_LR.wav',[I Q],fs);         % store result
% audiowrite('FMRadio_IQ_250kHz_L.wav',[I Q],fs);          % store result
  
  if(ifigs==1)
     m=1:5000; N = length(x); f = fs/N*(0:N-1);
     figure, plot(m,I(m),m,Q(m)); title('FM Signal Re/Im = IQ'); xlabel('n'); grid; pause;
     figure; plot(f,20*log10(abs(fft(x)))); xlabel('f (Hz)'); title('FFT of FM Signal'); grid; pause
     figure, spectrogram(x, 512, 512-64, 512, fs); title('STFT of FM Signal'); pause
  end
