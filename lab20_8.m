% lab20_8: IQk -> IQ0(n) -> IQn
% This is the complete, runnable script for Exercise 20.8.
% It includes the rectangular filter option and all necessary plotting.

% Set to 1 to run this program standalone (as per Exercise 20.8)
if(1) 
    clear all; close all;
    
    % =====================================================================
    % Parameters to Modify for Exercise 20.8
    % =====================================================================
    
    % Requirement 6: Change modulation type
    modtype = '4QAM'; % 2PAM, 4PAM, 8PAM, BPSK, QPSK, DQPSK, 8PSK, 4QAM, 16QAM
    
    % Requirement 7 & 10: Change filter type and roll-off
    psf_type = 'sqrt-rc'; % 'sqrt-rc', 'normal-rc', or 'rect' (for Req 7)
    r = 0.35;             % PS filter roll-off factor (Req 10)
    
    % Requirement 9: Change samples per symbol (K)
    K = 24; % samples per symbol (interpolation factor)
    
    % Requirement 5: Toggle between 0 and 1
    do_updown = 0; % 0/1 frequency up-down conversion
    % change from 0 to 1 does the 20.10 exercise
    
    % =====================================================================
    
    Nsymbols = 250; % number of symbols to transmit
    [IQcodes, Nstates, Nbits, R] = IQdef(modtype); % Get modulation definitions
    numbers = floor(Nstates*(rand(Nsymbols,1)-10*eps)); % symbol generation
    IQk = numbers2IQ(numbers, modtype, IQcodes); % Generate symbols
end

% --- Parameters from Listing 20.4 ---
do_figures = 1;

% Requirement 4: Set do_disturb = 0
do_disturb = 0; 

Ns = 8; % symbols per PS filter
fs = 240000; % sampling frequency in Hz
fcar = 50000; % carrier frequency in Hz

fcut = fs/(2*K); % PSF cut-off frequency
Npsf = Ns*K+1; Mpsf = (Npsf-1)/2; % PSF filter length and its half

% Add prefix/postfix to handle filter ramp-up/down
numbers = floor(Nstates*(rand(2*Ns,1)-10*eps)); % 2*Ns random carrier states
dummy = numbers2IQ(numbers, modtype, IQcodes);
IQdum = [dummy(1:Ns); IQk; dummy(Ns+1:end)]; % dummy prefix and postfix

IQ0 = zeros(1, length(IQdum)*K); % appending prefix & postfix
IQ0(1:K:end) = IQdum; % zero insertion

% --- Requirement 7: Add Rectangular Filter Option ---
fprintf('Generating pulse shaping filter (psf_type: %s)...\n', psf_type);
if(isequal(psf_type, 'rect'))
    % Requirement 7: Add the rectangular pulse shaping filter
    hpsf = [ones(1,K) zeros(1, (Ns-1)*K+1)] / K;
    % Make sure filter has the correct length
    hpsf(Npsf) = 0;
else
    % Original code for 'sqrt-rc' and 'normal-rc'
    if(isequal(psf_type, 'sqrt-rc'))
        rctype = 'sqrt';
    else
        rctype = 'normal';
    end
    hpsf = firrcos(Npsf-1, fcut, r, fs, 'rolloff', rctype);
end

% --- Transmitter Side ---
IQn = conv(IQ0, hpsf); 
IQn = IQn(Mpsf+1 : end-Mpsf); % pulse shaping in TX

if (do_figures)
    % --- Plot 1: Filter and TX Spectrum (Req 8) ---
    figure;
    
    % Plot filter impulse response
    subplot(2,2,1);
    plot(-(Mpsf):Mpsf, hpsf);
    title('PSF Impulse Response h(n)');
    xlabel('n'); grid on;
    
    % Plot filter frequency response
    subplot(2,2,2);
    plot_psd(hpsf, fs, 'Filter Frequency Response');
    
    % Plot TX Baseband Spectrum
    subplot(2,2,3);
    plot_psd(IQn, fs, 'TX Baseband Spectrum (IQn)');
    
    % Plot TX Constellation
    subplot(2,2,4);
    plot(real(IQn), imag(IQn), 'b.');
    hold on;
    plot(real(IQdum), imag(IQdum), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    title('TX Phasor Diagram (IQn)');
    xlabel('I(n)'); ylabel('Q(n)'); grid on; axis equal;
    legend('Interpolated Path', 'Symbol States');
    pause(0.1);
end


% --- Channel Simulation (Up/Down Conversion) ---
if (do_updown) % Frequency UP and DOWN, channel and disturbances are in between
    fprintf('Simulating Frequency UP/DOWN conversion...\n');
    % Frequency UP conversion in TX (quadrature modulator)
    N = length(IQn); n = 0:N-1; % signal length, sample indexes
    
    % This is the real, passband signal
    y = real(IQn).*cos(2*pi*fcar/fs*n) - imag(IQn).*sin(2*pi*fcar/fs*n);

    % ... Comments about channel/DAC/ADC from original listing ...
    
    % Frequency DOWN conversion in RX
    % (We assume no disturbances as per do_disturb = 0)
    df = 0; dphi = 0; % carrier frequency and phase offsets
    
    % This is the complex baseband signal at the receiver
    IQnn = 2*y.*exp(-j*(2*pi*(fcar/fs + df/fs)*n + dphi));
    
    if(do_disturb)
       % This block is skipped as per Exercise 20.8
    end
    
else
    % *** LOGIC FIX ***
    % If not doing up/down conversion, just pass the TX signal 
    % directly to the receiver's matched filter.
    % The original listing was missing this!
    fprintf('Bypassing Frequency UP/DOWN. Sending TX signal to RX filter.\n');
    IQnn = IQn;
end

% --- Receiver Side ---
% Apply matched filter (the same PSF)
IQnn = conv(IQnn, hpsf); 
IQnn = IQnn(Mpsf+1 : end-Mpsf); % pulse shaping in RX

% Synchronization and Sampling (ideal)
% We remove the prefix/postfix symbols and sample at the right time
IQnn_synced = IQnn(Ns*K+1 : end-Ns*K); % Remove ramp-up/down
IQkk = IQnn_synced(1:K:end); % Sample at symbol rate

if (do_figures)
    % --- Plot 2: Receiver Eye/Phasor Diagrams (Req 10) ---
    figure;
    
    % Plot RX Eye Diagram for I
    subplot(2,2,1);
    plot_eye(real(IQnn_synced), K, 'RX Eye Diagram - I(n)');
    
    % Plot RX Eye Diagram for Q
    subplot(2,2,2);
    plot_eye(imag(IQnn_synced), K, 'RX Eye Diagram - Q(n)');
    
    % Plot RX Phasor Diagram (paths)
    subplot(2,2,3);
    plot(real(IQnn_synced), imag(IQnn_synced), 'b.');
    title('RX Phasor Diagram (Paths)');
    xlabel('I(n)'); ylabel('Q(n)'); grid on; axis equal;
    
    % Plot RX Constellation (sampled points)
    subplot(2,2,4);
    plot(real(IQkk), imag(IQkk), 'r.', 'MarkerSize', 10);
    title('RX Constellation (Sampled Symbols)');
    xlabel('I(k)'); ylabel('Q(k)'); grid on; axis equal;
    
end

fprintf('Simulation complete. Compare plots to observe exercise requirements.\n');


% =====================================================================
% Helper Functions (from book archive / standard practice)
% =====================================================================

function [IQcodes, Nstates, Nbits, R] = IQdef(modtype)
    % Defines the constellation points (IQ codes) for modulations
    Nstates = 0; Nbits = 0; R = 1; IQcodes = [];
    if(isequal(modtype, '2PAM'))
        Nstates = 2; Nbits = 1; IQcodes = [-1, 1];
    elseif(isequal(modtype, '4PAM'))
        Nstates = 4; Nbits = 2; IQcodes = [-3, -1, 1, 3];
    elseif(isequal(modtype, '8PAM'))
        Nstates = 8; Nbits = 3; IQcodes = [-7, -5, -3, -1, 1, 3, 5, 7];
    elseif(isequal(modtype, 'BPSK'))
        Nstates = 2; Nbits = 1; IQcodes = [exp(j*0), exp(j*pi)];
    elseif(isequal(modtype, 'QPSK'))
        Nstates = 4; Nbits = 2; R = 1/sqrt(2); 
        IQcodes = [R+j*R, -R+j*R, -R-j*R, R-j*R]; % R*exp(j*pi/4*[1 3 5 7])
    elseif(isequal(modtype, 'DQPSK'))
        Nstates = 4; Nbits = 2; R = 1/sqrt(2);
        IQcodes = [R+j*R, -R+j*R, -R-j*R, R-j*R]; % Differential
    elseif(isequal(modtype, '8PSK'))
        Nstates = 8; Nbits = 3; 
        IQcodes = exp(j*pi/4*(0:7));
    elseif(isequal(modtype, '4QAM'))
        Nstates = 4; Nbits = 2; R = 1/sqrt(2);
        IQcodes = [R+j*R, -R+j*R, -R-j*R, R-j*R];
    elseif(isequal(modtype, '16QAM'))
        Nstates = 16; Nbits = 4; R = 1/sqrt(10);
        [I,Q] = meshgrid([-3, -1, 1, 3], [-3, -1, 1, 3]);
        IQcodes = R*(I + 1j*Q);
        IQcodes = IQcodes(:).'; % Gray coding would be better, but use book's
    else
        error('Unknown modulation type');
    end
end

function IQk = numbers2IQ(numbers, modtype, IQstates)
    % State numbers to IQ values
    if(isequal(modtype, 'DQPSK'))
        IQk(1) = 1; % Initial state
        for k = 1:length(numbers)
            IQk(k+1) = IQk(k) * IQstates(numbers(k)+1); % Differential coding
        end
    else
        IQk = IQstates(numbers + 1);
    end
    IQk = IQk(:); % Ensure column vector
end

function plot_psd(signal, fs, plot_title)
    % Helper function to plot Power Spectral Density
    
    % --- FIX START ---
    % The original hard-coded window 'hanning(1024)' fails if the 
    % signal is shorter than 1024 samples (like the hpsf filter).
    % This fix adjusts the window size to be no larger than the signal itself.
    
    L = length(signal);
    win_size = min(L, 1024); % Use a window size <= signal length
    n_overlap = floor(win_size / 2); % 50% overlap
    n_fft = 1024; % Keep 1024 for consistent FFT bin resolution
    
    [Pxx, F] = pwelch(signal, hanning(win_size), n_overlap, n_fft, fs, 'centered');
    % --- FIX END ---
    
    plot(F/1000, 10*log10(Pxx));
    title(plot_title);
    xlabel('Frequency (kHz)');
    ylabel('Power/Frequency (dB/Hz)');
    grid on;
end

function plot_eye(signal, K, plot_title)
    % Helper function to plot an eye diagram
    % K = samples per symbol
    len = floor(length(signal) / (2*K)) * (2*K);
    sig_matrix = reshape(signal(1:len), 2*K, []);
    
    % Plot 2 symbol durations
    t = (0:2*K-1) / K; 
    plot(t, sig_matrix, 'b');
    title(plot_title);
    xlabel('Symbols');
    ylabel('Amplitude');
    grid on;
end