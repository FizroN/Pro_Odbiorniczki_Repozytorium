% pluto_rx.m
% Receives and decodes a signal using the ADALM-PLUTO SDR.
% Based on lab21_ex_receiver.m

clear all; close all;

%% 1. Parameters (MUST MATCH THE TRANSMITTER)
% --- SDR Parameters ---
f_carrier = 915e6;         % 915 MHz (Must match TX)
pluto_ip = 'ip:192.168.2.1'; % !!! CHANGE THIS to your Pluto's IP

% --- Waveform Parameters ---
modtype = 'QPSK';      % Must match TX
Rs = 57e3;             % Symbol Rate (Must match TX)
K = 8;                 % Samples per Symbol (Must match TX)
Ns = 8;
fs = Rs * K;           % Baseband Sample Rate (Must match TX)
r = 0.35;              % Roll-off factor
rctype = 'sqrt';
Ndata = 1000;          % Number of data symbols (Must match TX)

% --- Receiver Control Flags (from lab21) ---
% We enable all synchronization and correction
do_synchro = 2;    % 2 = use more robust xD(n), sD(n) correlation
do_cfoequ = 2;     % 2 = use polyfit for robust CFO estimation
do_chanequ = 1;    % 1 = use simple one-tap (gain/phase) equalizer

%% 2. Generate Local Header Reference (from lab21)
% The receiver MUST know what the header looks like
% so it can search for it.
Npsf = K*Ns+1; Mpsf = (Npsf-1)/2;
[IQcodes, Nstates, Nbits, R ] = IQdef( modtype );
[numHead, Nhead ] = modtype2header( modtype );
IQkHead = numbers2IQ( numHead, modtype, IQcodes );
numData = floor( Nstates*(rand(Ndata,1)-10*eps) ); % For error check

% --- Generate the "ideal" pulse-shaped header ---
% This is 's(n)' from the original script
[IQnHead, hpsf_tx ] = IQ2psf( IQkHead, K, Ns, r, rctype );
s = IQnHead; 

% This is the "ideal" RRC filter for the TX side.
% The RX filter 'hpsf' will be the same.
hpsf = hpsf_tx;

%% 3. Set up ADALM-PLUTO Receiver
disp('Setting up Pluto RX...');

% --- Calculate how many samples to grab ---
% We need at least one full frame + header to find the signal.
% Let's grab 3 frames' worth to be safe.
frame_len = (Nhead + Ndata) * K;
samples_to_grab = frame_len * 3; 

try
    rx = sdrrx('Pluto', ...
        'RadioID', pluto_ip, ...
        'CenterFrequency', f_carrier, ...
        'BasebandSampleRate', fs, ...
        'SamplesPerFrame', samples_to_grab, ...
        'OutputDataType', 'double', ...
        'GainSource', 'AGC Slow Attack'); % AGC is helpful!

    %% 4. Capture Signal
    disp('Capturing signal...');
    % Flush the SDR buffer by capturing a few times
    for i = 1:5
        rx_signal_raw = rx();
    end
    
    % This is our captured data!
    % It will be a column vector
    rx_signal_raw = rx();
    disp('Signal captured!');
    release(rx); % We can release the hardware now

    if (isempty(rx_signal_raw))
        error('No signal captured. Check antennas and TX.');
    end

    % --- The received signal 'x(n)' ---
    % Make it a row vector to match the original script's format
    x = rx_signal_raw.';

    %% 5. Run the Receiver DSP Chain (from lab21)
    
    % --- Low-pass pulse shaping filter (MATCHED FILTER) ---
    disp('Applying matched filter...');
    x = conv( x, hpsf ); % 'hpsf' is from step 2

    % --- Differential signals for robust correlation ---
    sD = s(2:end) .* conj(s(1:end-1));  % sD(n) 
    xD = x(2:end) .* conj( x(1:end-1)); % signal xD(n)

    % --- Cross correlation function ---
    disp('Correlating to find header...');
    if( do_synchro == 1)   % less robust
       Cxs = conv( x-mean(x), conj(s(end:-1:1)-mean(s)) );
       [peak_val, nmax ] =  max( abs( Cxs ) );
       n1st = nmax - length(s) + (Npsf-1)/2 + 1;
    else % (do_synchro == 2) more robust
       CxDsD = conv( xD-mean(xD), conj(sD(end:-1:1)-mean(sD)) );
       [peak_val, nmax ] =  max( abs( CxDsD ) );
       n1st = nmax - length(sD) + (Npsf-1)/2 + 2;
    end
    
    fprintf('Synchronization found! Peak = %.2f, Index = %d\n', peak_val, n1st);
    
    % Plot correlation
    figure;
    if (do_synchro == 1)
        plot(abs(Cxs), '.-'); title('Cross-Correlation |Cxs(k)|');
    else
        plot(abs(CxDsD), '.-'); title('Cross-Correlation |CxDsD(k)|');
    end
    grid on; xlabel('k');
    hold on; plot(nmax, peak_val, 'rs', 'MarkerSize', 10);
    
    % --- Check if sync is valid ---
    if (n1st <= 0 || n1st > length(x) - frame_len)
        error('Synchronization failed! Found header at invalid index. Try running again.');
    end

    % --- Use header for CFO and Phase Estimation ---
    disp('Estimating frequency and phase offsets...');
    nsynch = n1st : K : n1st+(Nhead-1)*K; % Symbol indices of RX header
    nhead = Mpsf+1 : K : Mpsf+1+(Nhead-1)*K; % Symbol indices of local header
    
    % Compare received header to ideal header
    work = x( nsynch ) .* conj( s( nhead ) );
    
    if( do_cfoequ == 1 )      % simple estimator     
      df = mean( angle( conj( work(1:end-1) ) .* work(2:end) ) );
      df = df / (2*pi*K); dph = 0;
    end
    if( do_cfoequ == 2 )      % polyfit (robust)
      phi0 = angle(work(round(Nhead/2))); work=work.*exp(-j*phi0);
      ang = unwrap( angle(work) ); nn = 0 : K : (Nhead-1)*K; 
      temp = polyfit( nn, ang, 1); df = temp(1)/(2*pi); dph = temp(2)+phi0;
      
      figure; plot( nn,ang,'b.-',nn,temp(2)+temp(1)*nn,'r-'); grid;
      title('Phase vs. Time (Polyfit for CFO)'); xlabel('Sample'); ylabel('Phase (rad)');
      legend('Measured Phase', 'Fit Line');
    end
    
    fprintf('CFO Estimate = %.5f (normalized), Phase Estimate = %.2f rad\n', df, dph);

    % --- Apply CFO and Phase Correction ---
    disp('Applying correction...');
    % Cut the signal to start from the beginning of the header
    IQn_corrected = x( n1st : end ) .* exp(-j*(2*pi*df*(0:length(x)-n1st)+dph));
    
    % --- Channel Equalization (One-tap) ---
    % Use the (now corrected) header to find channel gain
    IQkHeadEstim = IQn_corrected( 1 : K : 1+(Nhead-1)*K );
    
    if( do_chanequ == 1 ) 
      % Find the average complex gain (attenuation + phase)
      gains = IQkHeadEstim .* conj(IQkHead) ./ abs(IQkHead).^2;
      gain = mean(gains);
      IQn_corrected = IQn_corrected / gain; % Apply correction
      fprintf('Channel Gain Estimate = %.2f, Phase = %.2f rad\n', abs(gain), angle(gain));
    end
    
    %% 6. Decode and Check Errors
    disp('Decoding data...');
    % Get the symbols from the corrected signal
    % We need to grab (Nhead + Ndata) symbols
    num_symbols_to_decode = Nhead + Ndata;
    if (length(IQn_corrected) < num_symbols_to_decode * K)
        warning('Not enough signal captured to decode full frame. Errors will be high.');
        num_symbols_to_decode = floor(length(IQn_corrected) / K);
    end
    
    ns_decode = 1 : K : (num_symbols_to_decode * K);
    IQk_received = IQn_corrected(ns_decode);
    
    % --- Final Constellation Plot ---
    figure;
    plot( real(IQk_received), imag(IQk_received), 'b.' );
    hold on;
    % Plot the ideal constellation points for reference
    plot(real(IQcodes), imag(IQcodes), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    grid on; title('RX Decoded Constellation (Corrected)');
    axis equal;
    
    % --- Convert symbols to numbers ---
    received_numbers = IQ2numbers( IQk_received, modtype );
    
    % --- Compare TX vs RX ---
    % We need to compare the *data* part, which is after the header
    rx_data = received_numbers(Nhead+1 : end).';
    tx_data = numData(1 : length(rx_data)); % Trim to match received length
    
    errors = sum( rx_data ~= tx_data );
    ber = errors / length(tx_data);
    
    fprintf('\n--- RESULTS ---\n');
    fprintf('Errors: %d / %d\n', errors, length(tx_data));
    fprintf('Bit Error Rate (BER): %.4f %%\n', ber * 100);

catch ME
    disp('Error setting up or running Pluto RX:');
    disp(ME.message);
    if (exist('rx', 'var'))
        release(rx);
    end
end
