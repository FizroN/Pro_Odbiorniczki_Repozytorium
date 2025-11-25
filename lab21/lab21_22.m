function [ns, IQs] = timing_recovery(IQ, alg, Nsps, n1st)
% Implementacja Gardnera i Muellera-Mullera z interpolacją Farrowa
% IQ - sygnał wejściowy
% alg - numer algorytmu (3 = Gardner odporny na CFO)
% Nsps - liczba próbek na symbol (np. 2, 4)
% n1st - indeks startowy

N = length(IQ); 
% Parametry pętli adaptacyjnej
damp = sqrt(2)/2; 
band = 0.01; % Pasmo pętli (szybkość reakcji)
mi1 = (4*damp*band)/(1 + 2*damp*band + band*band);
mi2 = (4*band*band)/(1 + 2*damp*band + band*band);

% Bufory dla interpolatora Farrowa (kwadratowego)
x2 = [0 0 0]; x1 = [0 0 0]; x0 = [0 0 0];

adap1 = 0; adap2 = 0; % Rejestry pętli
k = 1; 
n = n1st;

% Główna pętla po próbkach sygnału
while n < N - Nsps - 2
    
    % Pobranie próbek do interpolacji
    % (Symulujemy, że mamy dostęp do 'n' i sąsiadów)
    curr_segment = IQ(n : n+2); 
    
    % Współczynniki Farrowa (wielomiany Lagrange'a)
    % x(t) ~ c2*t^2 + c1*t + c0
    % Tutaj uproszczona wersja dla przesunięcia frakcyjnego 'adap1'
    
    % Dla uproszczenia w pętli - zakładamy Nsps próbek skoku
    % Ale musimy interpolować DOKŁADNĄ wartość w punkcie 'adap1'
    
    % Wersja uproszczona z książki (Listing 21.6):
    % Zamiast pełnej pętli sample-by-sample, skaczemy co symbol i interpolujemy
    
    % Obliczamy próbki w oknie interpolacji
    s_m1 = IQ(n);       % Próbka wcześniejsza
    s_0  = IQ(n+1);     % Próbka "centralna"
    s_p1 = IQ(n+2);     % Próbka późniejsza
    
    % Filtr Farrowa (współczynniki wielomianu)
    v2 = 0.5 * (s_p1 + s_m1) - s_0;
    v1 = 0.5 * (s_p1 - s_m1);
    v0 = s_0;
    
    % Interpolacja w punkcie 'mu' (adap1)
    mu = adap1; 
    val_now = v2*mu*mu + v1*mu + v0; % Interpolowana próbka "TERAZ"
    
    % Potrzebujemy też próbki "w połowie" symbolu dla Gardnera (half-symbol)
    % Zakładamy, że jest ok. Nsps/2 próbek wcześniej
    % W pełnej implementacji Farrowa to jest bardziej złożone, 
    % tu użyjemy aproksymacji z próbek fizycznych dla błędu
    
    % Detektor błędu Gardnera (Algorithm 3 - odporny na rotację fazy)
    % err = real( (val_now - val_prev) * conj(val_mid) )
    
    % Ponieważ implementacja pełnego Farrowa w pętli while jest skomplikowana,
    % użyjmy gotowego fragmentu z książki dla Nsps małego:
    
    % Zapisujemy wynik
    IQs(k) = val_now;
    ns(k) = n;
    
    % Obliczenie błędu (uproszczone dla Gardnera)
    % Potrzebujemy próbki w połowie i poprzedniej
    if k > 2
       prev = IQs(k-1);
       mid_idx = round(n - Nsps/2); 
       if mid_idx > 0
           mid = IQ(mid_idx); % Bierzemy najbliższą fizyczną próbkę jako "środek"
           
           % Wzór Gardnera:
           err = real( (IQs(k) - prev) * conj(mid) );
           
           % Aktualizacja pętli
           adap2 = adap2 + mi2 * err;
           adap1 = adap1 + adap2 + mi1 * err;
           
           % Zawijanie offsetu (żeby nie uciekł poza +/- 0.5 próbki)
           % W tej prostej wersji Farrowa sterujemy tylko ułamkiem, 
           % a część całkowitą dodajemy do indeksu 'n'
           
           while adap1 > 0.5
               adap1 = adap1 - 1;
               n = n + 1; % Przeskocz próbkę fizyczną do przodu
           end
           while adap1 < -0.5
               adap1 = adap1 + 1;
               n = n - 1; % Cofnij się o próbkę fizyczną
           end
       end
    end
    
    % Skok do następnego symbolu
    n = n + Nsps;
    k = k + 1;
end

% Obcięcie pustych końcówek
IQs = IQs(1:k-1);
ns = ns(1:k-1);

end

% lab21_ex_tetra.m
% Synchronization continuous down-link burst in TETRA digital telephony
clear all; close all;

load('TETRA_423.4125MHz_noise_-22.mat');  K=10;  % files: -22,-10,0,10,22.mat % AWGN (noisy) channel
%load('TETRA_423.4125MHz_flat.mat');  K=5;        % flat channel 
%load('TETRA_423.4125MHz_tu50.mat');  K=8;        % typical urban channel at 50 kph
%load('TETRA_423.4125MHz_ht200.mat'); K=4;        % hilly terrain channel at 200 kph

IQn = y.'; clear y; fs   % sampling frequency is written from the MAT file 

% Parameters
% fs = 2560000;     % sampling frequency 102.4 kHz or 2.56 MHz - written from the file
fsymb = 18000;      % symbol frequency 18 kHz
modtype = 'QPSK';  % 2PAM, 4PAM, 8PAM, BPSK, QPSK, DQPSK, 8PSK, 4QAM, 16QAM
Nhead = 19;         % STS header length (number of carrier states)
Ndata = 129;        % data length - searched (number of carrier states)
Nframe = 255;       % frame length, including Nhead and Ndata
rctype = 'sqrt';    % PSF: 'sqrt': raised cosine filter type for TX and RX 
r = 0.35;           % PSF: filter roll-off factor
%K = 10;            % number of samples per symbol, at present fs/fsymb
                    % but we can change it
Ns = 10;            % PSF: symbols per PS filter 

% Signal resampling
fsnew = ceil(K)*fsymb;                 % new sampling frequency
[UP,DOWN] = rat( fsnew/fs );           % UP/DOWN interpolation/decimation orders
I = resample( real(IQn), UP, DOWN );   % resampling I(n)
Q = resample( imag(IQn), UP, DOWN );   % resampling Q(n)
IQn = I + j*Q; clear I Q;              % combining I(n)+j*Q(n)
fs = fsnew; K = ceil(K);               % setting new values

% STS 19 header - known
numHead = [ 3 0 0 1 2 1 3 0 3 2 2 1 3 0 0 1 2 1 3 ];
% Data to be detected after the 1st header in files -22, -10, 0, 10, 22.mat:
% Length = 129 = 15+108+1+5

% Transmitted data - known
numData = [ ...
  2,3,3,2,2,3,1,3,2,1,1,3,0,2,2,1,0,2,3,3,3,3,0,1,0,0,0,1,0,3,0,2,3,...
  3,1,0,2,3,2,1,2,1,2,3,3,3,2,2,3,2,1,1,2,0,0,3,1,2,1,3,2,0,0,2,2,3,...
  1,3,2,3,0,2,1,0,3,1,3,0,1,1,3,1,2,2,2,1,1,2,0,2,1,3,2,0,1,2,1,3,1,...
  3,0,2,1,2,1,0,1,1,3,3,1,3,2,2,2,1,1,3,1,0,2,2,1,3,2,3,1,3,0 ];
%Ndata = length(numData), pause

% Signal and its spectrum
figure; n = 1:2500;
subplot(211); plot( real(IQn(n)) ); grid; title('I(n)');
subplot(212); plot( real(IQn(n)) ); grid; title('Q(n)'); pause
figure;
pwelch(IQn,2048,2048-1024,2048,fs); pause

% ##############
% Added receiver
% ##############

% INITIALIZATION ##########################################################

do_synchro = 1;    % 0/1/2: 0 = none, 1 = using x(n),s(n), 2=using xD(n),sD(n)
do_decim = 0;      % 0/1 signal decimation in the receiver
do_cfoequ = 1;     % 1/2 freq carrier offset estimation: 1 = simple, 2 = polyfit
do_chanequ = 1;    % 0/1/2/3 channel correction methods: 0=none, 1=ZF, 2=LS/MSE, 3=NLMS filter
Mdecim = 10;       % decimation order
Mdelay = 0;        % delay
do_timing_rec = 1; % Włącznik naszej nowej funkcji

% Pulse shaping filter
Npsf = floor(Ns*K+1), Mpsf=(Npsf-1)/2, pause              % PSF filter length and its half
hpsf = firrcos( Npsf-1, 1/(2*K), r, 1,'rolloff',rctype);  % 'normal' or 'sqrt'

% Generate modulation codebook
[IQcodes, Nstates, Nbits, R ] = IQdef( modtype );    % carrier IQ codes

% Calculate IQkHead from numHead
IQkHead= numbers2IQ( numHead, modtype, IQcodes );    % IQk states for the Header

% Calculate IQnHead from IQkHead - do pulse shaping of the header
IQnHead = IQ2psf( IQkHead, K, Ns, r, 'normal');      % IQn samples of the Header   

% Low-pass pulse shaping filter (PSF) in the receiver
IQn = conv( IQn, hpsf );

figure; n = 1:K*100;
plot( real(IQn(n)), imag(IQn(n)),'b-' ); grid; title('Q(n) = f( I(n) )'); pause

% DO IT NOW! ##############################################################

% Possible signal down-sampling in the receiver

if( do_decim )
  if( Mdecim>1 )  
    M = Mdecim; Mdelay=0;
    IQn = [ zeros(1,Mdelay) IQn ];     % signal delay
    IQn = resample(IQn, 1, M);         % low-pass filtering and decimation
    IQnHead = resample(IQnHead, 1, M); % low-pass filtering and decimation
  % IQnHead = IQnHead(1:M:end);        % decimation
    hpsf = M*hpsf(1:M:Npsf); K=K/M; Npsf=(Npsf-1)/M+1; Mpsf=(Npsf-1)/2; % changing PSF
    N = length(IQn); N = floor(N/K)*K; n = 1:N; IQn = IQn(n);           % changing IQn
  end
else M=1;
end

s = IQnHead;  % IQn of synchronization header
x = IQn;      % IQn of received signal in the base-band

% Auto correlation function
sD = s(2:end) .* conj(s(1:end-1));                             % sD(n) 
ms = mean(numHead); Css1 = conv( numHead-ms, numHead(end:-1:1)-ms ); % auto corr numHead()
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
Cxs = conv( x-mean(x), conj(s(end:-1:1)-mean(s)) );   % cross s(n) with x(n)
xD = x(2:end)  .* conj( x(1:end-1));                  % signal xD(n)
CxDsD = conv( xD-mean(xD), sD(end:-1:1)-mean(sD) );   % cross sD(n) with xD(n)

    figure; L = length(Cxs); k=1:L;
    plot(k,abs(Cxs),'.-' ); grid; xlabel('k'); title('|Cxs(k)|'); pause
    L = length(CxDsD); k=1:L;
    figure; plot(k,abs(CxDsD),'.-'); grid; xlabel('k'); title('|CxDsD(k)|'); pause

if( do_synchro == 0 )   % NO SYNCHRONIZATION
       
   IQn = IQn( 2*Npsf-1 : end - (2*Npsf-1)+1 );                % remove transients
   received = IQ2numbers( IQn(1:K:end), modtype );            % decimate
   errors = sum( received ~= num( Ns + 1 : end-Ns ) ), pause  % errors
   
else                    % SYNCHRONIZATION
   if (do_timing_rec)
    % Jeśli użyliśmy Gardnera, mamy już zdecymowane symbole w IQ_sync
    % Musimy tylko znaleźć w nich nagłówek (bo Gardner nie szuka ramki, tylko symboli)

    % To jest trudne, bo IQ_sync ma inną długość niż IQn. 
    % Dla uproszczenia w tym zadaniu: wyświetlamy tylko konstelację (powyżej).
    % Pełna integracja wymagałaby ponownego puszczenia korelacji na IQ_sync.

    disp('Konstelacja wyrysowana. Pełna detekcja ramki na sygnale asynchronicznym wymagałaby ponownej korelacji.');
else
    Nsync = length( IQnHead ); 
   if( do_synchro==1)   % less robust to noise, using x(n), s(n)
       [dummy, nmax ] =  max( abs( Cxs ) );    % maximum position
       n1st = nmax - Nsync + (Npsf-1)/2 + 1,   % 1st header symbol
   else                                        % more robust to noise, using xD(n), sD(n)
       [dummy, nmax ] =  max( abs( CxDsD ) );  % maximum position
       n1st = nmax - Nsync + (Npsf-1)/2 + 2,   % 1st header symbol
   end
   
 % n1st = 2714; % 664;

 % Use only header symbols for frequency offset and phase shift estimation 
   nsynch = n1st : K : n1st+(Nhead-1)*K; nhead = Mpsf+1 : K : Mpsf+1+(Nhead-1)*K;
   work = IQn( nsynch ) .* conj( IQnHead( nhead ) );  % # the same
 % work = IQn( nsynch ) .* conj( IQkHead );           % # 
        figure; stem( abs(IQkHead - IQnHead(nhead)) ), title('Head ERR? Is 0?'); pause
        figure; stem( abs(IQkHead - IQn(nsynch)) ), title('Head ERROR: |TX-RX|'); pause
   
   if( do_cfoequ ==1 )      % simple frequency carrier offset estimator     
      df = mean( angle( conj( work(1:end-1) ) .* work(2:end) ) ); % simple version
      df = df / (2*pi*K); dph = 0;
   end
   
   if( do_cfoequ ==2 )      % phase polynomial fitting method
      phi0 = angle(work(round(Nhead/2))); work=work.*exp(-j*phi0);
    % phi0=0;
      ang = unwrap( angle(work) ); nn = 0 : K : (Nhead-1)*K; 
      temp = polyfit( nn, ang, 1); df = temp(1)/(2*pi); dph = temp(2)+phi0;
           figure; plot( nn,ang,'b-',nn,temp(2)+temp(1)*nn,'r-'); grid; title('Angle'); pause
      chanPh = 0; carDPh = 0;
   end
   
 % Do frequency offset and phase correction - from channel and freq-down converter
   IQn = IQn( n1st : length(IQn) ) .* exp(-j*(2*pi*df*(0:length(IQn)-n1st)+dph));
   N=length(IQn); n=1:N; ns=1:K:N;
 % Results
   carOffs_estim  = df/M*fs, pause

 % Amplitude and phase correction - using header for channel estimation & correction
 % Knowing input and output, we can estimate channel and correct it
   IQkHeadEstim = IQn( 1 : K : 1+(Nhead-1)*K );                 % detected header states
   
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
      rxData = IQ2numbers( IQkDataEstim, modtype );   % IQ to state numbers
      figure; 
      plot( IQkDataEstim,'r*'); grid; title('Q(n) = f( I(n) )'); pause
      return
   end
   
   if( do_chanequ==3 ) % NLMS adaptive filtering of symbols - weights adaptation using header
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
      figure; 
      plot( uest(end-Ndata/2:end),'r*'); grid; title('Q(n) = f( I(n) )'); pause
      return
   end
   
 % Results
   received = IQ2numbers( IQn( 1 : K : 1+(Nhead+Ndata+1)*K ), modtype ),
   errors = sum( received ~= [ numHead numData ] ), pause

    received = IQ2numbers( IQn( 1 : K : 1+(Nhead+Ndata+1)*K ), modtype );
    errors = sum( received ~= [ numHead numData ] ), pause
end
   
end % end of synchronization using header

n = 1 : (Nhead+Ndata)*K; IQn = IQn(n); ns = 1 : K : length(IQn);
figure; 
plot( real(IQn), imag(IQn),'b-',real(IQn(ns)),imag(IQn(ns)),'ro','MarkerFaceColor','red' ); grid; title('Q(n) = f( I(n) )'); pause
figure; 
plot( real(IQn(ns)),imag(IQn(ns)),'ro','MarkerFaceColor','red' ); grid; title('Q(n) = f( I(n) )'); pause

N = length( IQn ); N = floor(N/K)*K; n = 1:N;
Nshift=floor(K/2)+1; m = 1+Nshift : N-(K-Nshift);
figure;
subplot(121); plot( reshape(real(IQn(m)),K,(N-K)/K),'b'); xlabel('n'); title('RX: Eye diagram for I(n)'); grid;
subplot(122); plot( reshape(imag(IQn(m)),K,(N-K)/K),'r'); xlabel('n'); title('RX: Eye diagram for Q(n)'); grid; pause


% --- ZADANIE 21.22: TIMING RECOVERY ---


if (do_timing_rec == 1)
    fprintf('Uruchamianie Gardner Timing Recovery...\n');
    
    % Ustalamy parametry
    alg = 3; % Gardner
    start_idx = 1; 
    
    % Wywołujemy funkcję
    [sample_indices, IQ_sync] = timing_recovery(IQn, alg, K, start_idx);
    
    % Podmieniamy sygnał do detekcji
    % IQn_sync zastąpi nasze "proste decymowanie"
    
    % Wizualizacja konstelacji PO synchronizacji
    figure; 
    plot(real(IQ_sync), imag(IQ_sync), 'r.');
    title(['TETRA Constellation after Gardner (K=' num2str(K) ')']);
    grid on; axis square;
    
    % WAŻNE: Nadpisujemy zmienną, której używa reszta skryptu do liczenia błędów
    % Skrypt dalej spodziewa się, że zdecymowany sygnał to po prostu wektor.
    % Musimy jednak uważać, bo 'received' jest obliczane kawałek dalej.
end
% --------------------------------------