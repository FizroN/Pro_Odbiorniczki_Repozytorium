% lab21_ex_receiver.m
% Testing synchronization procedures in single carrier receiver
clear all; close all;

modtype = '4QAM';  % 2PAM, 4PAM, 8PAM, BPSK, QPSK, DQPSK, 8PSK, 4QAM, 16QAM
rctype = 'sqrt';   % 'normal, 'sqrt': raised cosine filter type for TX and RX 
Ndata = 250;       % number of carrier states (symbols) to be transmitted
K = 24; Ns = 8;    % PSF: samples per symbol, symbols per pulse shaping filter 
r = 0.35;          % PSF: filter roll-off factor
fs = 240000;       % sampling frequency in Hz: 1, 250e+3, 240e+3
fcar = 40000;      % carrier frequency in Hz

do_updownchan = 0; % 0/1 frequency up-down conversion plus channel simulation
do_disturb = 1;    % 0/1 addition of disturbances in the baseband
do_synchro = 1;    % 0/1/2: 0 = none, 1 = using x(n),s(n), 2=using xD(n),sD(n)
do_decim = 0;      % 0/1 signal decimation in the receiver
do_cfoequ = 1;     % 0/1/2 freq carrier offset estimation: 0=none, 1=simple, 2=polyfit
do_chanequ = 1;    % 0/1/2/3 channel correction methods: 0=none, 1=ZF, 2=LS/MSE, 3=NLMS filter

%chan = [ 0.5, -0.25, 0.1, -0.1 ]; % channel impulse response in baseband symbol-spaced
chan = [ 1 ];                      % perfect channel
SNR=160;  chanG=1; chanPh=0;  carDF=0.0000; carDPh=0;  ADCdt=0;           % No disturb
%SNR=40;  chanG=0.25; chanPh=pi/7; carDF=0.0002; carDPh=pi/11; ADCdt=0.5; % Disturb
Mdecim = 1;        % decimation order: 1, 2, 3, 4, 6, 8, 12, 24
Mdelay = 0;        % decimation delay: in samples before decimation  

% Old transmission part - as before
Npsf = K*Ns+1; Mpsf = (Npsf-1)/2;
[IQcodes, Nstates, Nbits, R ] = IQdef( modtype );       % take carrier IQ codes
% IQk of Header
[numHead, Nhead ] = modtype2header( modtype );          % take header IQ numbers
IQkHead= numbers2IQ( numHead, modtype, IQcodes );       % calculate IQ states
% IQk of Data
numData = floor( Nstates*(rand(Ndata,1)-10*eps) );      % generate random IQ numbers
IQkData = numbers2IQ( numData, modtype, IQcodes );      % calculate IQ states
% Numbers ALL, IQk ALL
num = [ numData' numHead' numData' numHead' numData' ]; % ALL transmitted IQ numbers
IQk = [ IQkData  IQkHead  IQkData  IQkHead  IQkData ];  % ALL transmitted IQ states
% IQn of Header only (pulse shaping)
IQnHead = IQ2psf( IQkHead, K, Ns, r, 'normal');         % IQn of Header   
% IQn of everything (pulse shaping)
[IQn, hpsf ] = IQ2psf( IQk, K, Ns, r, rctype );         % IQn of ALL

do_timing_rec = 1; 

if (do_timing_rec)
    % Wybór algorytmu:
    % 1 - Gardner (wrażliwy na offset freq)
    % 2 - Gardner (na składowej Q)
    % 3 - Gardner (odporny na offset freq)
    % 4 - Gardner (wersja sprzężona)
    % 5 - Mueller-Muller (wymaga 1 próbki na symbol, reszta 2+)
    alg = 1; 
    
    % Mdecim musi być ustawione tak, aby K/Mdecim (próbki na symbol) >= 2 dla Gardnera
    % lub >= 1 dla Muellera-Mullera.
    % Dla przykładu, jeśli K=24, a chcemy 2 próbki na symbol, Mdecim = 12.
    
    % Pętla po różnych opóźnieniach startowych (dla testu)
    for Mdelay = 1 : 1 : K
        
        % Przygotowanie sygnału wejściowego (z decymacją wstępną)
        % Wybieramy co Mdecim-tą próbkę
        IQnn = IQn(1+Mdelay : Mdecim : end);
        
        nist = 1; % Indeks startowy
        
        % WYWOŁANIE FUNKCJI TIMING RECOVERY
        % Funkcja ta (opisana w Listingu 21.6) musi być dostępna w katalogu roboczym!
        % [nowe_indeksy, wartosci_symboli] = funkcja(sygnał, algorytm, Nsps, start)
        [ns, IQs] = timing_recovery(IQnn, alg, K/Mdecim, nist);
        
        % Tu normalnie nastąpiłaby analiza błędów dla 'IQs'
        % W ramach ćwiczenia program wygeneruje wykresy wewnątrz funkcji timing_recovery
    end
end

    N = length( IQn ); n = Npsf : N-Npsf+1; ns = Npsf : K : N-Npsf+1;
    figure; plot( real(IQn(n)), imag(IQn(n)), real(IQn(ns)), imag(IQn(ns)),'ro','MarkerFaceColor','red'); grid; title('TX: Q(n) = f( I(n) )'); pause
    n = n + floor(K/2)+1;
    figure
    subplot(121); plot( reshape(real(IQn(n)),K,length(n)/K),'b'); xlabel('n'); title('TX: Eye diagram for I(n)'); grid;
    subplot(122); plot( reshape(imag(IQn(n)),K,length(n)/K),'r'); xlabel('n'); title('TX: Eye diagram for Q(n)'); grid; pause

if(0) % 0/1 Testing symbol timing recovery methods    
    alg = 1; % 1/2/3/4/5 timing recovery algorithm
    for Mdelay = 1 : 1 : K
        Mdelay
        IQnn = IQn( 1+Mdelay : Mdecim : end ); n1st = 1;
        [ ns, IQs ] = timing_recovery( IQnn, alg, K/Mdecim, n1st );
    end
end

% Optional frequency UP and DOWN conversion in TX plus channel simulation  
if( do_updownchan )
   if( length( chan ) > 1 )         % when more than one channel/filter weight
      chan = resample(chan,K,1);    % upsampling channel impulse response
      figure; plot(chan); title('h(n)'); pause
   end 
   df = 0; dphi = 0;  % CFO added in the base-band OR df = carDF*fs; dphi = carDPh;
   IQn = IQupchandown( IQn, fcar, fs, chan, df, dphi );
end

% Addition of disturbances in the base-band
if( do_disturb )
    IQn = IQdisturb( IQn, SNR, chanG, chanPh, carDF, carDPh, ADCdt, Npsf );
    
        N = length( IQn ); n = Npsf : N-Npsf+1; ns = Npsf : K : N-Npsf+1;
        figure; plot( real(IQn(n)), imag(IQn(n)), real(IQn(ns)), imag(IQn(ns)),'ro','MarkerFaceColor','red'); grid; title('RX: Q(n) = f( I(n) )'); pause
        n = n + floor(K/2)+1;
        figure
        subplot(121); plot( reshape(real(IQn(n)),K,length(n)/K),'b'); xlabel('n'); title('RX: Eye diagram for I(n)'); grid;
        subplot(122); plot( reshape(imag(IQn(n)),K,length(n)/K),'r'); xlabel('n'); title('RX: Eye diagram for Q(n)'); grid; pause
end

% Possible signal down-sampling in the receiver

if( do_decim )
  if( M > 1 )  
    M = Mdecim;                        % copy initial program parameter
    IQn = [ zeros(1,Mdelay) IQn ];     % signal delay
    IQn = resample(IQn, 1, M);         % low-pass filtering and decimation
    IQnHead = resample(IQnHead, 1, M); % low-pass filtering and decimation
  % IQnHead = IQnHead(1:M:end);        % only decimation
    hpsf = M*hpsf(1:M:Npsf); K=K/M; Npsf=(Npsf-1)/M+1; Mpsf=(Npsf-1)/2;
    N = length(IQn); N = floor(N/K)*K; n = 1:N; IQn = IQn(n);
  end
else M=1;
end

% Low-pass pulse shaping filter in the receiver
IQn = conv( IQn, hpsf );

s = IQnHead;  % IQn of synchronization header
x = IQn;      % IQn of received signal in the base-band

% Auto correlation function
sD = s(2:end) .* conj(s(1:end-1));                             % sD(n) 
ms = mean(numHead); Css1 = conv( numHead-ms, numHead(end:-1:1)-ms ); % auto corr numHead 
ms = mean(s);       Css2 = conv( s-ms, conj(s(end:-1:1)-ms) ); % auto corr s(n)
ms = mean(sD);      Css3 = conv( sD-ms, sD(end:-1:1)-ms );     % auto corr sD(n)
    
    figure; stem( numHead ); grid; xlabel('n'); title('s(n)'); pause
    figure; L=length(Css1); L=(L-1)/2; k=-L:L;
    stem( k, Css1 ); grid; xlabel('k'); title('Css1(k)'); pause
    figure; L=length(Css2); L=(L-1)/2; k=-L:L;
    plot( k, abs(Css2),'.-' ); grid; xlabel('k'); title('|Css2(k)|'); pause
    figure; L=length(Css3); L=(L-1)/2; k=-L:L;;
    plot( k, abs(Css3),'.-' ); grid; xlabel('k'); title('|Css3(k)|'); pause

% Cross correlation function
Cxs = conv( x-mean(x), conj(s(end:-1:1)-mean(s)) );   % cross corr s(n) with x(n)
xD = x(2:end)  .* conj( x(1:end-1));                  % signal xD(n)
CxDsD = conv( xD-mean(xD), sD(end:-1:1)-mean(sD) );   % cross corr sD(n) with xD(n)

    figure; L = length(Cxs); k=1:L;
    plot(k,abs(Cxs),'.-' ); grid; xlabel('k'); title('|Cxs(k)|'); pause
    L = length(CxDsD); k=1:L;
    figure; plot(k,abs(CxDsD),'.-'); grid; xlabel('k'); title('|CxDsD(k)|'); pause

if( do_synchro == 0 )   % NO SYNCHRONIZATION
       
   IQn = IQn( 2*Npsf-1 : end - (2*Npsf-1)+1 );                % remove transients
   received = IQ2numbers( IQn(1:K:end), modtype );            % decimate
   errors = sum( received ~= num( Ns + 1 : end-Ns ) ), pause  % errors
   
else                    % SYNCHRONIZATION
    
   Nsync = length( IQnHead ); 
   if( do_synchro==1)   % less robust to noise, using x(n), s(n)
       [dummy, nmax ] =  max( abs( Cxs ) );    % maximum position
       n1st = nmax - Nsync + (Npsf-1)/2 + 1,   % 1st header symbol
   else                 % more robust to noise, using xD(n), sD(n)
       [dummy, nmax ] =  max( abs( CxDsD ) );  % maximum position
       n1st = nmax - Nsync + (Npsf-1)/2 + 2,   % 1st header symbol
   end

 % Use only header symbols for frequency offset and phase shift estimation 
   nsynch = n1st : K : n1st+(Nhead-1)*K; nhead = Mpsf+1 : K : Mpsf+1+(Nhead-1)*K;
   work = IQn( nsynch ) .* conj( IQnHead( nhead ) );  % # the same
 % work = IQn( nsynch ) .* conj( IQkHead );           % # 
        figure; stem( abs(IQkHead - IQnHead(nhead)) ), title('Head ERR? Is 0?'); pause
        figure; stem( abs(IQkHead - IQn(nsynch)) ), title('Head ERROR: |TX-RX|'); pause
   if( do_cfoequ == 0 ) df=0; dph = 0; end
   if( do_cfoequ == 1 )      % simple frequency carrier offset estimator     
      df = mean( angle( conj( work(1:end-1) ) .* work(2:end) ) ); % simple version
      df = df / (2*pi*K); dph = 0;
   end
   if( do_cfoequ == 2 )      % phase polynomial fitting method
      phi0 = angle(work(round(Nhead/2))); work=work.*exp(-j*phi0);
    % phi0=0;
      ang = unwrap( angle(work) ); nn = 0 : K : (Nhead-1)*K; 
      temp = polyfit( nn, ang, 1); df = temp(1)/(2*pi); dph = temp(2)+phi0;
           figure; plot( nn,ang,'b-',nn,temp(2)+temp(1)*nn,'r-'); grid;
           title('Angle'); pause
      if( carDF == 0) allPhase_estim = dph, expected = (chanPh+carDPh), end
      chanPh = 0; carDPh = 0;
   end
 % Do frequency offset and phase correction - from channel and freq-down converter
   IQn = IQn( n1st : length(IQn) ) .* exp(-j*(2*pi*df*(0:length(IQn)-n1st)+dph));
   N=length(IQn); n=1:N; ns=1:K:N;
 % Results
   carOffs_estim  = df/M,  expected = carDF, pause

 % Amplitude and phase correction - using header for channel estimation & correction
 % Knowing input and output, we can estimate a channel and correct it
   IQkHeadEstim = IQn( 1 : K : 1+(Nhead-1)*K );             % detected header states
   if( do_chanequ==0 ) gain=1; dph = 0; end                 % no corrector
   if( do_chanequ==1 ) % one-tap corrector 
      gains = IQkHeadEstim .* conj(IQkHead) ./ abs(IQkHead).^2; % compared with known
      gain = mean(real(gains)) + j*mean(imag(gains));           % mean channel "gain"
      IQn = IQn / gain;                                         % signal correction
   end
   if( do_chanequ==2 ) % solving linear input-output equation in LS (MSE) sense
      M = 9;                           % channel length (in number of symbols)
      L = Nhead;                       % number of input symbols
      v = IQkHeadEstim;                % values of output symbols
      V = toeplitz(v(M:L),v(M:-1:1));  % matrix with output signal 
      u = IQkHead(M:L).';              % input symbols
      g = (V\u),                       % channel inverse filter
      v = IQn( 1+(Nhead-1-(M-2))*K : K : 1 + (Nhead+Ndata-1)*K); % data to correct
      L = length(v);                   % number of input data
      V = toeplitz(v(M:L),v(M:-1:1));  % matrix with input data
      IQkDataEstim = V*g;              % data correction
      [ IQkHeadEstim(1:10).' IQkHead(1:10).' ], pause     % compare IQ header 
      [ IQkDataEstim(1:10)   IQkData(1:10).' ], pause     % compare IQ data
      rxData = IQ2numbers( IQkDataEstim, modtype );       % IQ to state numbers
      [ rxData(1:15).' numData(1:15) ],                   % numbers in and out
      errors = sum( rxData.' - numData ), pause           % error
      figure; 
      plot( IQkDataEstim,'r*'); grid; title('Q(n) = f( I(n) )'); pause
      return
   end  
   if( do_chanequ==3 ) % NLMS adaptive filter of symbols - weights adaptation using header
      u = [ IQkHead IQkData ];                              % sent
      v = IQn( 1 : K : 1+(Nhead+Ndata-1)*K );               % received
      bv = zeros(1,M); g = zeros(1,M); mi=0.1; ghist = [];  % initialization
      for n = 1 : length(u)                                 % filter loop
          bv = [ v(n) bv(1:M-1) ];                          % input to buffer
          uest(n) = sum( bv .* conj(g) );                   % estimated value
          err(n) = u(n) - uest(n);                          % error
          g = g + mi * conj(err(n)) * bv / (bv*bv');        % filter weights update
          ghist = [ ghist g(1) ];                           % history
      end                                                   %
      figure; plot( abs(ghist) ); title('ghist(n)'); pause  % figure
      rx = IQ2numbers( uest, modtype ).';                   % received state numbers
      tx = [ numHead; numData ];                            % transmitted state numbers 
      errors = sum( rx(end-Ndata/2:end) - tx(end-Ndata/2:end) ), pause  % error
      figure; 
      plot( uest(end-Ndata/2:end),'r*'); grid; title('Q(n) = f( I(n) )'); pause
      return
   end    
 % Results
   if( carDF == 0) allPhase_estim = angle(gain), expected = chanPh+carDPh, end
   chanGain_estim = abs(gain), expected = chanG, pause
 % Errors
   received = IQ2numbers( IQn( 1 : K : (Nhead+Ndata)*K ), modtype );
 % errors = sum( received' ~= [ numHead; numData ] ), pause
   errors = sum( received(Nhead+1:Nhead+Ndata)' ~= [ numData ] ), pause
end % end of synchronization using the header

% Final figures

n = 1 : (Nhead+Ndata)*K; IQn = IQn(n); ns = 1 : K : length(IQn);
figure; 
plot( real(IQn), imag(IQn),'b-',real(IQn(ns)),imag(IQn(ns)),'ro','MarkerFaceColor','red' ); grid; title('Q(n) = f( I(n) )'); pause

N = length( IQn ); N = floor(N/K)*K; n = 1:N;
Nshift=floor(K/2)+1; m = 1+Nshift : N-(K-Nshift);
figure;
subplot(121); plot( reshape(real(IQn(m)),K,(N-K)/K),'b'); xlabel('n'); title('RX: Eye diagram for I(n)'); grid;
subplot(122); plot( reshape(imag(IQn(m)),K,(N-K)/K),'r'); xlabel('n'); title('RX: Eye diagram for Q(n)'); grid; pause
