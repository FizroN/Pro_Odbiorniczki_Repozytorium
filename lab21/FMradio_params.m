% Ten plik musiałem trochę zedytować, żeby 21.11 zadziałał mi...
% fmradio_params.m 
% FM Radio - initialization of FM radio parameters and weights of filters
% In Octave: >>pkg load signal[ENTER]

% PARAMETER INITIALIZATION, FILTER DESIGN   
  fs = 250000;         % sampling frequency of one FM radio station
  fpilot = 19000;      % frequency of the pilot, 19000 Hz
  fsymb = fpilot/16;   % frequency of RDS symbols 1187.5 Hz, 19000/16
  fstereo = 2*fpilot;  % frequency of the L-R signal carrier, 38000 Hz
  frds = 3*fpilot;     % frequency of the RDS signal carrier, 57000 Hz
  faudio = 25000;      % frequency of an audio signal (assumed)
  L = 500;             % length of used FIR filters
  Ks = 6;              % number of symbols in the PSF filter
  dt = 1/fs;           % sampling ratio

% Pre-emphasis and de-emphasis filtes for frequency faudio
  f1 = 2120; tau1 = 1/(2*pi*f1); w1 = tan(1/(2*faudio*tau1));
  b_de = [w1/(1+w1), w1/(1+w1)]; a_de = [1, (w1-1)/(w1+1)];
  b_pre = a_de; a_pre = b_de;

% Pulse shaping filter (PSF) for RDS symbols
  Tsymb = 1/fsymb;           % time duration of one RDS symbol   
  Nsymb = round(fs/fsymb);   % number of samples per one RDS symbol
  Npsf = Ks*Nsymb;           % number of samples of the PSF filter
  if(rem(Npsf,2)==1) Npsf=Npsf+1; end
  df = fs/Npsf; f = 0 : df : 2/Tsymb; Nf=length(f);        % spectrum 
  H = zeros(1, Npsf); H(1:Nf) = cos(pi*f*Tsymb/4); H(end:-1:end-(Nf-2)) = H(2:Nf);
  hpsf1 = fftshift(ifft(H)); hpsf1=hpsf1/max(hpsf1);       % imp. response #1
  hpsf2 = firrcos( Npsf, fsymb, 1.0, fs,'rolloff','sqrt'); % imp. response #2
  hpsf2 = hpsf2/max(hpsf2); hpsf2 = hpsf2(1:end-1);
  n1=1:length(hpsf1); n2=1:length(hpsf2);
% plot(n1,hpsf1,'ro-',n2,hpsf2,'bx-'); grid; title('hpsf1(n) (R) and hpsf2(n) (B)'); pause
  hpsf = hpsf2;
  phasePSF = angle( exp(-j*2*pi*fsymb/fs*[0:Npsf-1]) * hpsf' ); % phase shift

% Low-pass (LP) filter for recovery of the L+R signal
  hLPaudio = fir1(L,(faudio/2)/(fs/2),kaiser(L+1,7));
% Narrow band-pass (BP) filter for separation of the pilot signal (around 19 kHz)
  fcentr = fpilot;  df1 = 1000; df2 = 2000;
  ff = [ 0 fcentr-df2  fcentr-df1 fcentr+df1 fcentr+df2 fs/2 ]/(fs/2);
  fa = [ 0 0.01        1          1          0.01       0 ];
  hBP19 = firls(L,ff,fa);  % Matlab, Octave; firpm() only for Matlab
% --- POPRAWKA DLA FILTRU 38 kHz ---
% Używamy fir1, bo jest stabilniejszy i daje poprawne wzmocnienie (0 dB)
% Pasmo przepustowe: 38 kHz +/- 10 kHz (szeroko, żeby złapać całe audio L-R)
Wn_stereo = [28000 48000] / (fs/2); % Normalizacja do Nyquista
hBP38 = fir1(L, Wn_stereo, 'bandpass');
% ----------------------------------
% hBP38 = firls(L,ff,fa);  % Matlab, Octave; firpm() only for Matlab  
% Narrow band-pass filter for separation of the RDS signal (around 57 kHz)
  fcentr = frds; df1 = 2000; df2 = 4000;
  ff = [ 0 fcentr-df2  fcentr-df1 fcentr+df1 fcentr+df2 fs/2 ]/(fs/2);
  fa = [ 0 0.01        1          1          0.01       0    ];
  hBP57 = firls(L,ff,fa);  % Matlab, Octave; firpm() only for Matlab 
% Narrow band-pass filter for separation of the 2*fstereo component (around 76 kHz)
  fcentr = 2*fstereo; df1 = 0.8*100; df2 = 100;
  ff = [ 0 fcentr-df2  fcentr-df1 fcentr+df1 fcentr+df2 fs/2 ]/(fs/2);
  fa = [ 0 0.01        1          1          0.01       0    ];
  hBP76 = firls(L,ff,fa);  % Matlab, Octave; firpm() only for Matlab