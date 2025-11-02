% pluto_loopback_test.m
%
% Performs a full-duplex (simultaneous TX/RX) loopback test on a 
% single ADALM-PLUTO SDR.
%
% *** HARDWARE SETUP ***
% 1. Connect one antenna to the PLUTO TX port.
% 2. Connect a second antenna to the PLUTO RX port.
% 3. Place the antennas near each other (e.g., 30 cm / 1 ft apart).
%
clear all; close all;

%% 1. Parameters (Combined TX and RX)
% --- SDR Parameters ---
f_carrier = 915e6;              % 915 MHz (Unlicensed ISM band)
pluto_ip = 'ip:192.168.2.1';     % Default IP for USB connection
tx_gain = -25;                  % *** FIX 1: Increased gain from -30 to -20 ***
rx_gain_mode = 'AGC Slow Attack'; % Let receiver auto-adjust its gain

% --- Waveform Parameters (MUST BE CONSISTENT) ---
modtype = 'QPSK';      
Rs = 57e3;                  % Symbol Rate (57 kS/s)
K = 8;                      % Samples per Symbol
Ns = 8;                     % Symbols per pulse shaping filter
fs = Rs * K;                % Baseband Sample Rate (456 kS/s)
r = 0.35;                   % Roll-off factor
rctype = 'sqrt';            % Root-Raised Cosine
Ndata = 1000;               % Number of data symbols

% --- Receiver Control Flags ---
do_synchro = 2;    % 2 = use robust xD(n), sD(n) correlation
do_cfoequ = 2;     % 2 = use polyfit for robust CFO estimation
do_chanequ = 1;    % 1 = use simple one-tap (gain/phase) equalizer

%% 2. Generate TX Frame and Header Reference
disp('Generating waveform and receiver references...');
Npsf = K*Ns+1; Mpsf = (Npsf-1)/2;
[IQcodes, Nstates, Nbits, R ] = IQdef( modtype );

% --- Create Header and Data Symbols ---
[numHead, Nhead ] = modtype2header( modtype );
IQkHead = numbers2IQ( numHead, modtype, IQcodes );
numData = floor( Nstates*(rand(Ndata,1)-10*eps) );
IQkData = numbers2IQ( numData, modtype, IQcodes );

% --- Combine into a single frame ---
IQk = [ IQkHead  IQkData ];
num = [ numHead' numData' ]; % Full sequence of transmitted numbers (for error check)

% --- Pulse Shape the TX Signal ---
[IQn, hpsf ] = IQ2psf( IQk, K, Ns, r, rctype );

% --- Prepare TX signal for SDR (column vector, scaled) ---
tx_signal = IQn / max(abs(IQn));
tx_signal = tx_signal.'; % Must be a column vector

% --- Generate the "ideal" header 's' for the receiver's correlator ---
[s, ~] = IQ2psf( IQkHead, K, Ns, r, rctype );

frame_len_samples = length(tx_signal);
fprintf('Waveform generated. Frame length: %d samples.\n', frame_len_samples);

%% 3. Setup SDR (TX and RX)
disp('Setting up Pluto TX and RX objects...');

try
    % --- Configure Transmitter ---
    tx = sdrtx('Pluto', ...
        'RadioID', pluto_ip, ...
        'CenterFrequency', f_carrier, ...
        'BasebandSampleRate', fs, ...
        'Gain', tx_gain); % Use the stronger gain

    % --- Configure Receiver ---
    % *** FIX 2: Capture 1 full second (fs samples) to guarantee catching
    %     the short transmission. This is a much larger 'net'. ***
    samples_to_grab = fs; 
    
    rx = sdrrx('Pluto', ...
        'RadioID', pluto_ip, ...
        'CenterFrequency', f_carrier, ...
        'BasebandSampleRate', fs, ...
        'SamplesPerFrame', samples_to_grab, ...
        'OutputDataType', 'double', ...
        'GainSource', rx_gain_mode);
    
    fprintf('Receiver buffer size set to %d samples (1.0 second).\n', samples_to_grab);

    %% 4. Transmit and Receive (The Loopback)
    disp('*** Transmitting and Receiving Simultaneously... ***');
    
    % FLUSH THE RX BUFFER FIRST
    disp('Flushing RX buffer...');
    for i=1:5
        rx(); 
    end 
    
    % Transmit the signal
    disp('Transmitting signal...');
    tx(tx_signal); 
    
    % Grab the 1-second-long buffer
    disp('Grabbing 1-second capture buffer...');
    rx_signal_raw = rx();
    
    disp('*** Capture complete! ***');
    
    %% 5. Release Hardware
    release(tx);
    release(rx);
    disp('Hardware released.');

    %% 6. Run the Receiver DSP Chain
    
    % --- Use the captured signal as 'x' ---
    x = rx_signal_raw.'; % Make it a row vector
    
    % --- Apply Matched Filter ---
    disp('Applying matched filter...');
    x = conv( x, hpsf ); 

    % --- Differential signals for robust correlation ---
    sD = s(2:end) .* conj(s(1:end-1));  
    xD = x(2:end) .* conj( x(1:end-1)); 

    % --- Cross correlation ---
    disp('Correlating to find header in 1-second buffer...');
    CxDsD = conv( xD-mean(xD), conj(sD(end:-1:1)-mean(sD)) );
    [peak_val, nmax ] =  max( abs( CxDsD ) );
    n1st = nmax - length(sD) + (Npsf-1)/2 + 2;
    
    % *** NEW CHECK: See if the peak is actually a peak ***
    if (peak_val < 1e-3) % A threshold to check if we found *anything*
       error('Correlation peak is still 0.00! Signal is too weak. Try moving antennas closer or increasing tx_gain to -10.'); 
    end
    
    fprintf('Synchronization found! Peak = %.2f, Index = %d\n', peak_val, n1st);
    
    % Plot correlation
    figure;
    plot(abs(CxDsD), '.-'); 
    title('Cross-Correlation |CxDsD(k)|'); grid on; xlabel('k');
    hold on; plot(nmax, peak_val, 'rs', 'MarkerSize', 10);
    
    % --- Check if sync is valid ---
    if (n1st <= 0 || n1st + frame_len_samples > length(x))
        error('Synchronization failed! Found header at invalid index.');
    end

    % --- CFO and Phase Estimation ---
    disp('Estimating frequency and phase offsets...');
    nsynch = n1st : K : n1st+(Nhead-1)*K; 
    nhead = Mpsf+1 : K : Mpsf+1+(Nhead-1)*K;
    work = x( nsynch ) .* conj( s( nhead ) );
    
    % Polyfit (robust)
    phi0 = angle(work(round(Nhead/2))); work=work.*exp(-j*phi0);
    ang = unwrap( angle(work) ); nn = 0 : K : (Nhead-1)*K; 
    temp = polyfit( nn, ang, 1); df = temp(1)/(2*pi); dph = temp(2)+phi0;
      
    figure; 
    plot( nn,ang,'b.-',nn,temp(2)+temp(1)*nn,'r-'); grid;
    title('Phase vs. Time (Polyfit for CFO)'); xlabel('Sample'); ylabel('Phase (rad)');
    legend('Measured Phase', 'Fit Line');
    fprintf('CFO Estimate = %.5f (normalized), Phase Estimate = %.2f rad\n', df, dph);

    % --- Apply CFO and Phase Correction ---
    disp('Applying correction...');
    IQn_corrected = x( n1st : end ) .* exp(-j*(2*pi*df*(0:length(x)-n1st)+dph));
    
    % --- Channel Equalization (One-tap) ---
    IQkHeadEstim = IQn_corrected( 1 : K : 1+(Nhead-1)*K );
    gains = IQkHeadEstim .* conj(IQkHead) ./ abs(IQkHead).^2;
    gain = mean(gains);
    IQn_corrected = IQn_corrected / gain; % Apply correction
    fprintf('Channel Gain Estimate = %.2f, Phase = %.2f rad\n', abs(gain), angle(gain));
    
    %% 7. Decode and Check Errors
    disp('Decoding data...');
    num_symbols_to_decode = Nhead + Ndata;
    if (length(IQn_corrected) < num_symbols_to_decode * K)
        warning('Not enough signal captured to decode full frame.');
        num_symbols_to_decode = floor(length(IQn_corrected) / K);
    end
    
    ns_decode = 1 : K : (num_symbols_to_decode * K);
    IQk_received = IQn_corrected(ns_decode);
    
    % --- Final Constellation Plot ---
    figure;
    plot( real(IQk_received), imag(IQk_received), 'b.' ); % Received
    hold on;
    plot(real(IQcodes), imag(IQcodes), 'ro', 'MarkerSize', 10, 'LineWidth', 2); % Ideal
    grid on; title('RX Decoded Constellation (Corrected)');
    axis equal;
    
    % --- Convert symbols to numbers ---
    received_numbers = IQ2numbers( IQk_received, modtype );
    
    % --- Compare TX vs RX ---
    rx_data = received_numbers(Nhead+1 : end).';
    tx_data = numData(1 : length(rx_data)); 
    
    errors = sum( rx_data ~= tx_data );
    ber = errors / length(tx_data);
    
    fprintf('\n--- RESULTS ---\n');
    fprintf('Errors: %d / %d\n', errors, length(tx_data));
    fprintf('Bit Error Rate (BER): %.4f %%\n', ber * 100);

catch ME
    disp('*** ERROR OCCURRED ***');
    disp(ME.message);
    % Make sure hardware is released on error
    if (exist('tx', 'var'))
        release(tx);
    end
    if (exist('rx', 'var'))
        release(rx);
    end
end

