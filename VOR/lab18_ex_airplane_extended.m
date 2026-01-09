% lab18_ex_airplane.m
  clear all; close all;
  m=128; cm_plasma=plasma(m); cm = plasma;   % color maps for gray printing

% Read a recorded IQ signal - two VOR avionics signals
 % FileName = 'SDRSharp_Airplane_112500kHz_IQ.wav'; T=5; demod=3; fc = 2.9285e+5;
 FileName = 'SDRSharp_Airplane_112849kHz_IQ.wav'; T=5; demod=3; fc = -5.5353e+4;

% #########################################################################
% START - From program lab16_ex_IQ_DFT.m ##################################
  inf = audioinfo(FileName), pause                % what is "inside" 
  fs = inf.SampleRate;                            % sampling rate
  if(T==0) [x,fs] = audioread(FileName);          % read the whole signal
  else     [x,fs] = audioread(FileName,[1,T*fs]); % read only T seconds 
  end                                             %
  whos, pause                                     % what is in the memory
  Nx = length(x),                                 % signal length

% Reconstruct the complex-value IQ data, if necessary add Q=0
  [dummy,M] = size(x);
  if(M==2) x = x(:,1) - j*x(:,2); else x = x(1:Nx,1) + j*zeros(Nx,1); end 
      nd = 1:2500;
      figure(1); plot(nd,real(x(nd)),'bo-',nd,imag(x(nd)),'r*--'); xlabel('n'); grid; 
      title('I(n) = (o) BLUE/solid    |    Q(n)= (*) RED/dashed'); pause

% Parameters - central signal sample, lengths of FFT and STFT
  Nc = floor( Nx/2 ); Nfft = min(2^17,2*Nc); Nstft = 512; 

% Power Spectral Density (PSD) of the signal
  n = Nc-Nfft/2+1 : Nc+Nfft/2;             % FFT length, samples for analysis
  df = fs/Nfft;                            % df - step in frequency
  f = df * (0 : 1 : Nfft-1);               % frequency axis [ 0, fs ]
  fshift = df * (-Nfft/2 : 1 : Nfft/2-1);  % frequency axis [ -fs/2, fs/2 ]
  w = kaiser(Nfft,10);                     % window function used
  X = fft( x(n) .* w );                    % DFT of windowed signal
  P = 2*X.*conj(X) / (fs*sum(w.^2));       % Power Spectral Density (dB/Hz)
  Pshift = fftshift( P );                  % circularly shifted PSD
  
  % Parameters for Short Time Fourier Transform (STFT) of the signal
  df = fs/Nstft; ff = df*(0:1:Nstft-1);  ffshift = df*(-Nstft/2:1:Nstft/2-1);

      figure(2)
      subplot(211); plot(f,10*log10(abs(P))); xlabel('f (HZ)'); ylabel('(dB/Hz)')
      axis tight; grid; title('PSD for frequencies [0-fs)');
      subplot(212); spectrogram(x(n),kaiser(Nstft,10),Nstft-Nstft/4,ff,fs);
      colormap(cm); title('Short-time Fourier transform of x(n)'); pause

      figure(3)
      subplot(211); plot(fshift,10*log10(abs(Pshift))); xlabel('f (HZ)'); ylabel('(dB/Hz)')
      axis tight; grid; title('PSD for frequencies [-fs/2, fs/2)');
      subplot(212); spectrogram(x(n),kaiser(Nstft,10),Nstft-Nstft/4,ffshift,fs);
      colormap(cm); title('Short-time Fourier transform of x(n)'); pause
      subplot(111);

 % STOP - From program lab16_ex_IQ_DFT.m ##################################
 % ########################################################################
 
 if(demod==3) % airplane azimuth decoding using VOR signal

    % --- SYSTEM CONFIGURATION ---
    % M: Filter order. Higher = sharper filtering.
    M = 501; M2=(M-1)/2;                
    % fam: Frequency bandwidth (25 kHz is standard for VOR avionics).
    fam = 25000; dt=1/fam;              
    f1 = fc-fam/2; f2 = fc+fam/2;       
    
    %%% EDUCATIONAL: FILTER DESIGN
    % We design a band-pass filter to isolate the VOR signal from the raw radio spectrum.
    h = remez(M-1,[0 (fam/2-df) (fam/2+df) fs/2]/(fs/2),[ 1 0.9 0.1 0 ],[1 1]); h=reshape(h,1,M); h=h.*exp(j*2*pi*fc/fs*(-M2:M2)); % Matlab (old), Octave 

    % --- VOR SIGNAL ISOLATION ---
    x = conv(x,h); x=x(M:end-M+1);      % Apply filter and compensate delay
    
    %%% EDUCATIONAL: AM DEMODULATION (ENVELOPE DETECTION)
    % VOR uses Amplitude Modulation (AM) for the 30Hz variable signal.
    % We calculate the Envelope using the Pythagorean theorem on I/Q components.
    x = sqrt( real(x).*real(x) + imag(x).*imag(x) ); 
    
    %%% EDUCATIONAL: DOWNSAMPLING (DECIMATION)
    % The raw signal is sampled at ~1 MHz. The audio content is < 15 kHz.
    % We reduce the sample rate to 25 kHz (fam) to:
    % 1. Speed up calculations.
    % 2. Convert radio signal to "audio" signal range.
    x = decimate(x, round(fs/fam));     
    x = x - mean(x);                    % Remove DC offset

    % !!! CAPTURE POINT 1: We save the AM demodulated signal.
    % This contains both the 30Hz tone and the 9960Hz subcarrier.
    x_AM_demod = x; 
    xc = x; % Make a backup copy for parallel processing
    
    % =========================================================================
    % PATH 1: EXTRACTING THE VARIABLE PHASE SIGNAL (30 Hz AM)
    % =========================================================================
    
    %%% EDUCATIONAL: LOW-PASS FILTERING
    % We need to extract the pure 30 Hz sine wave.
    % We use a Low-Pass filter with cutoff ~50 Hz to remove everything else.
    hLP30 = fir1(M-1, 50/(fam/2), 'low');   
    x = conv(x,hLP30); x = x(M2+1:end-M2);  
    x = x - mean(x);                        
    
    % Normalization using max(abs()) is robust against signal inversion.
    x_azim = x(2:end-1)/max(abs(x));        
    
    % =========================================================================
    % PATH 2: EXTRACTING THE REFERENCE PHASE SIGNAL (FM on 9960 Hz)
    % =========================================================================
    
    %%% EDUCATIONAL: SUB-CARRIER ISOLATION
    % The Reference 30Hz signal is hidden inside a 9960Hz subcarrier (FM modulated).
    % First, we isolate this 10kHz region using a Band-Pass filter.
    hBP10k = fir1(M-1, [9000, 11000]/(fam/2), 'bandpass'); 
    x = conv(xc, hBP10k); x=x(M2+1:end-M2); 
    
    % !!! CAPTURE POINT FOR EXERCISE 18.11 ("Straightening the wheel")
    % We capture the signal before and during phase processing to visualize "Unwrap".
    x_subcarrier = x;                   % The pure 9.96 kHz sine wave (modulated)
    x_analytic = hilbert(x_subcarrier); % Complex analytic signal
    phase_raw = angle(x_analytic);      % "Sawtooth" phase (-pi to pi)
    
    %%% EDUCATIONAL: FM DEMODULATION (PHASE UNWRAPPING)
    % This step corresponds to Exercise 18.11 in the book.
    % The raw phase jumps from -pi to pi. We must "unwrap" it to make it continuous.
    % Only then can we calculate the derivative (Frequency).
    x = unwrap(phase_raw);              % "Straightening the wheel"
    phase_unwrapped = x;                % Save for plotting
    
    % Calculate Instantaneous Frequency (Derivative of Phase)
    x = x(3:end)-x(1:end-2);            % 3-point numerical differentiation
    x = (1/(2*pi))*x/(2*dt);            % Scale to Hz
    
    x_freq_noisy = x;                   % Save noisy FM signal
    
    %%% EDUCATIONAL: SMOOTHING
    % The raw FM signal is noisy. We know the payload is a 30Hz sine,
    % so we apply the same 30Hz Low-Pass filter as before.
    x = conv(x,hLP30); x=x(M2+1:end-M2);    
    x = x - mean(x);                        
    x_ref = x/max(abs(x));                  % Normalize

    % =========================================================================
    % RESULTS & VISUALIZATION
    % =========================================================================
    
    %%% EDUCATIONAL: AZIMUTH CALCULATION
    % Azimuth is simply the phase difference between the two 30Hz signals.
    phi_inst = angle( hilbert(x_azim) .* conj(hilbert(x_ref)) );
    phi_estim = mean( phi_inst );
    
    fprintf('\n--- VOR DECODING RESULTS ---\n');
    fprintf('Phase Shift (Radians): %.4f\n', phi_estim);
    fprintf('Airplane Azimuth:      %.2f degrees\n', rad2deg(phi_estim));
    fprintf('----------------------------\n');

    % --- FIGURE 1: TECHNICAL DASHBOARD (Classic Figure 18.8) ---
    figure('Name', 'VOR Technical Analysis', 'Color', 'white');
    
    subplot(3,2,1);
    % We use FFT to show the spectrum. Note the peak at 9.96kHz!
    Nfft = 4096; f_axis = (-Nfft/2 : Nfft/2-1) * (fam/Nfft);
    plot(f_axis, 20*log10(abs(fftshift(fft(x_AM_demod, Nfft)))));
    xlim([-12000 12000]); grid on; title('PSD of Signal (Note 9.96kHz subcarrier)');
    xlabel('Frequency (Hz)'); ylabel('dB');

    subplot(3,2,2);
    % This shows raw FM demodulation. It looks like noise but hides the signal.
    t_axis_short = (0:length(x_freq_noisy)-1)*dt;
    plot(t_axis_short, x_freq_noisy); xlim([2 2.1]); grid on;
    title('Raw FM Output (Before Smoothing)'); xlabel('Time (s)');

    subplot(3,2,3);
    % Spectrogram proves that 9.96kHz tone is "waving" (FM modulated).
    spectrogram(x_AM_demod, 256, 200, 1024, fam, 'yaxis'); ylim([0 12]);
    title('Spectrogram (Wavy line = FM)');

    subplot(3,2,4);
    % THE RESULT: Comparing the two sine waves. Shift = Direction.
    t_axis_final = (0:length(x_azim)-1)*dt;
    idx = find(t_axis_final > 2 & t_axis_final < 2.2);
    plot(t_axis_final(idx), x_azim(idx), 'b', 'LineWidth', 1.5); hold on;
    plot(t_axis_final(idx), x_ref(idx), 'r--', 'LineWidth', 1.5);
    legend('Variable (AM)', 'Reference (FM)'); grid on;
    title('Phase Comparison (Shift = Azimuth)');

    subplot(3,2,6);
    % Stability of the calculation over time.
    plot(t_axis_final(idx), phi_inst(idx)); yline(phi_estim, 'g--', 'Mean');
    grid on; title('Calculated Phase Shift [rad]');

    % --- FIGURE 2: EDUCATIONAL DEEP DIVE (Based on Ex 18.11 & Pilot View) ---
    figure('Name', 'Educational Extras', 'Color', 'white');
    
    % 1. Compass Plot (Pilot's View)
    subplot(1,2,1);
    polarplot([0 phi_estim], [0 1], 'r-^', 'LineWidth', 3);
    title({'COCKPIT VIEW', ['Bearing: ' num2str(rad2deg(phi_estim), '%.1f') ' deg']});
    rlim([0 1]); rticks([]); % Clean look
    
    % 2. The "Unwrap" Mystery (Exercise 18.11)
    subplot(1,2,2);
    t_zoom = 100:300; % Zoom in on a few samples
    plot(t_zoom, phase_raw(t_zoom), 'r.-'); hold on;
    plot(t_zoom, phase_unwrapped(t_zoom) - phase_unwrapped(t_zoom(1)), 'b-');
    legend('Raw Phase (Sawtooth)', 'Unwrapped Phase (Continuous)');
    title('Exercise 18.11: Phase Unwrapping');
    grid on; xlabel('Sample'); ylabel('Phase (radians)');
    
    % --- AUDIO DEMO  ---
    % fprintf('Playing VOR Audio (25kHz sample rate)...\n');
    % x_audio = x_AM_demod / max(abs(x_AM_demod));
    % volume = 0.55;
    % sound(x_audio * volume, fam); % Play the demodulated AM signal
    
    pause;
end
