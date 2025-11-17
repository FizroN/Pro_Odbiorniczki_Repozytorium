% lab20_ex_text_tx.m
% This script is for Exercise 20.15: 64-QAM and 256-QAM.
% It sends a text message through the full transmission chain,
% allowing for observation of disturbances and eye diagrams.

clear all; close all;

% =====================================================================
% 1. Load Text Payload
% =====================================================================
fprintf('Loading text payload...\n');
text_in = ['This is a longer text message for testing transmission. ', ...
           'We will use 64-QAM and 256-QAM to see how many errors ', ...
           'appear when noise is added to the channel. 0123456789. ', ...
           'The quick brown fox jumps over the lazy dog.'];
payload_len_bytes = length(text_in);

fprintf('Loaded text message: "%s"\n', text_in);
fprintf('Total bytes: %d\n', payload_len_bytes);

% =====================================================================
% 2. Simulation Parameters
% =====================================================================
% --- Modulation Parameters (Change for 20.15) ---
modtype = '256QAM'; % 4QAM, 16QAM, 64QAM, 256QAM

% --- Filter Parameters ---
psf_type = 'sqrt-rc'; % 'sqrt-rc', 'normal-rc', or 'rect'
r = 0.35;             % PS filter roll-off factor
K = 24; % samples per symbol (interpolation factor)
Ns = 8; % symbols per PS filter

% --- Channel Parameters ---
do_updown = 1; % 0/1 frequency up-down conversion
do_disturb = 1; % 0/1 to enable disturbances

fs = 240000; % sampling frequency in Hz
fcar = 50000; % carrier frequency in Hz

% --- Disturbance Parameters (Exercise 20.12) ---
% Ideal values (no disturbance):
%SNR_dB = 160; chan_gain = 1; chan_phase = 0; carDf = 0; carDph = 0; ADCppm = 0;

% Disturbed values (uncomment to test "rainy days"):
SNR_dB = 40; chan_gain = 1; chan_phase = 0; carDf = 0; carDph = 0; ADCppm = 0;

fprintf('--- Using Disturbances ---\n');
fprintf('SNR: %d dB\n', SNR_dB);
fprintf('Modulation: %s\n', modtype);


% =====================================================================
% 3. Transmitter Side
% =====================================================================
fprintf('Transmitting...\n');
[IQcodes, Nstates, Nbits, R] = IQdef(modtype);

% --- Convert payload (text) into modulation symbols (numbers) ---
[numbers, ~] = text2numbers(text_in, Nbits);
fprintf('Packed %d bytes into %d symbols using %d-bits/symbol.\n', ...
        payload_len_bytes, length(numbers), Nbits);

% Generate IQ symbols
IQk = numbers2IQ(numbers, modtype, IQcodes);

% --- Pulse Shaping ---
fcut = fs/(2*K); 
Npsf = Ns*K+1; Mpsf = (Npsf-1)/2; 
N_prefix = Ns;
numbers_dummy = floor(Nstates*(rand(2*N_prefix,1)-10*eps));
dummy = numbers2IQ(numbers_dummy, modtype, IQcodes);
IQdum = [dummy(1:N_prefix); IQk; dummy(N_prefix+1:end)]; 

IQ0 = zeros(1, length(IQdum)*K); 
IQ0(1:K:end) = IQdum; 

if(isequal(psf_type, 'rect'))
    hpsf = [ones(1,K) zeros(1, (Ns-1)*K+1)] / K; hpsf(Npsf) = 0;
else
    rctype = 'sqrt';
    if isequal(psf_type, 'normal-rc'), rctype = 'normal'; end
    hpsf = firrcos(Npsf-1, fcut, r, fs, 'rolloff', rctype);
end

IQn = conv(IQ0, hpsf); 
IQn = IQn(Mpsf+1 : end-Mpsf); 

% =====================================================================
% 4. Channel Simulation
% =====================================================================

if (do_updown) 
    fprintf('Simulating Frequency UP/DOWN conversion...\n');
    N = length(IQn); n = 0:N-1;
    
    y = real(IQn).*cos(2*pi*fcar/fs*n) - imag(IQn).*sin(2*pi*fcar/fs*n);
    
    df = 0; dphi = 0; 
    IQnn = 2*y.*exp(-j*(2*pi*(fcar/fs + df/fs)*n + dphi));
    
    if(do_disturb)
        fprintf('Applying disturbances...\n');
        IQnn = IQnn * chan_gain * exp(j*chan_phase);
        s_power = mean(abs(IQnn).^2);
        if s_power == 0, s_power = 1; end 
        noise_power = s_power / (10^(SNR_dB/10));
        noise = sqrt(noise_power/2) * (randn(size(IQnn)) + 1j*randn(size(IQnn)));
        IQnn = IQnn + noise;
        carDf_total = carDf + (ADCppm * 1e-6 * fcar / fs) / (1 + ADCppm * 1e-6);
        IQnn = IQnn .* exp(j*(2*pi*carDf_total*n + carDph));
    end
else
    fprintf('Bypassing Frequency UP/DOWN.\n');
    IQnn = IQn;
end

% =====================================================================
% 5. Receiver Side
% =====================================================================
fprintf('Receiving...\n');
% Apply matched filter
IQnn = conv(IQnn, hpsf); 
IQnn = IQnn(Mpsf+1 : end-Mpsf); 

% Synchronization and Sampling
IQnn_synced = IQnn(N_prefix*K+1 : end-N_prefix*K); 
IQkk_received = IQnn_synced(1:K:end); 

% --- Plot Eye/Phasor Diagrams (as requested by 20.15) ---
figure;
subplot(1,2,1);
plot_eye(real(IQnn_synced), K, sprintf('RX Eye Diagram - I(n) - %s', modtype));
subplot(1,2,2);
plot(real(IQkk_received), imag(IQkk_received), '.', 'MarkerSize', 5);
title(sprintf('RX Constellation (SNR: %d dB)', SNR_dB));
xlabel('I(k)'); ylabel('Q(k)'); grid on; axis equal;


% --- Decode IQ symbols back to numbers ---
received_numbers = IQ2numbers(IQkk_received, modtype, IQcodes);

% --- Convert numbers back to text ---
text_out = numbers2text(received_numbers, Nbits);

fprintf('Received and decoded %d bytes.\n', length(text_out));

% =====================================================================
% 6. Verification (Display Text)
% =====================================================================

% --- Trim output to match input length for comparison ---
text_out_trimmed = text_out(1:min(length(text_out), length(text_in)));

% --- Check for errors ---
num_errors = sum(text_in ~= text_out_trimmed);
BER = num_errors * 8 / (length(text_in) * 8); % Bit Error Rate (approx)

fprintf('\n--- TRANSMISSION COMPLETE ---\n');
fprintf('Original text:  %s\n', text_in);
fprintf('Received text:  %s\n', text_out_trimmed);
fprintf('\n');
fprintf('Modulation: %s\n', modtype);
fprintf('SNR: %d dB\n', SNR_dB);
fprintf('Character Errors: %d / %d\n', num_errors, length(text_in));
fprintf('Approximate Bit Error Rate (BER): %f\n', BER);
fprintf('Done.\n');


% =====================================================================
% HELPER FUNCTIONS
% =====================================================================

function [numbers, bitsnum] = text2numbers(text, Nbits)
    % text to IQ state numbers conversion
    bitschar = dec2bin(double(text), 8)'; % text-array, letters in rows
    [rows, cols] = size(bitschar);
    N = rows * cols; % number of all bits
    work = reshape(bitschar, [1, N]); % bits in one row
    Nadd = Nbits - rem(N, Nbits);
    if Nadd == Nbits, Nadd = 0; end % Handle case where N is multiple of Nbits
    for k = 1:Nadd, work = [work, '0']; end % appending 0 bits at the end
    bitsnum = reshape(work, [Nbits, (N+Nadd)/Nbits])'; % bits of all states
    numbers = bin2dec(bitsnum); % state numbers: from 0 to 2^Nbits-1
end

function text = numbers2text(numbers, Nbits)
    % IQ state numbers to text conversion
    text_bits = dec2bin(numbers, Nbits);
    [rows, cols] = size(text_bits);
    text_stream = reshape(text_bits', [1, rows*cols]); % one big stream of chars '0' '1'
    N = length(text_stream); 
    N = N - rem(N, 8); % Find last full byte
    text_stream = text_stream(1:N); % remove appended bits
    text_bytes = reshape(text_stream, [8, N/8])'; % strings of bytes
    text = char(bin2dec(text_bytes))'; % conversion to text
end


function [IQcodes, Nstates, Nbits, R] = IQdef(modtype)
    Nstates = 0; Nbits = 0; R = 1; IQcodes = [];
    if(isequal(modtype, '2PAM'))
        Nstates = 2; Nbits = 1; IQcodes = [-1, 1];
    elseif(isequal(modtype, '4PAM'))
        Nstates = 4; Nbits = 2; IQcodes = [-3, -1, 3, 1]; % Gray
    elseif(isequal(modtype, '8PAM'))
        Nstates = 8; Nbits = 3; IQcodes = [-7, -5, -3, -1, 7, 5, 3, 1]; % Gray
    elseif(isequal(modtype, 'BPSK'))
        Nstates = 2; Nbits = 1; IQcodes = [1, -1];
    elseif(isequal(modtype, 'QPSK'))
        Nstates = 4; Nbits = 2; R = 1/sqrt(2); 
        IQcodes = [R+j*R, -R+j*R, -R-j*R, R-j*R]; 
    elseif(isequal(modtype, 'DQPSK'))
        Nstates = 4; Nbits = 2; R = 1;
        IQcodes = [exp(j*pi/4), exp(j*3*pi/4), exp(j*7*pi/4), exp(j*5*pi/4)];
    elseif(isequal(modtype, '8PSK'))
        Nstates = 8; Nbits = 3; 
        idx = [0, 1, 3, 2, 7, 6, 4, 5]; % Gray
        IQcodes(idx+1) = exp(j*pi/4*(0:7));
    elseif(isequal(modtype, '4QAM'))
        Nstates = 4; Nbits = 2; R = 1/sqrt(2);
        IQcodes = [R*(-1-j*1), R*(-1+j*1), R*(1+j*1), R*(1-j*1)]; % Gray
    elseif(isequal(modtype, '16QAM'))
        Nstates = 16; Nbits = 4; R = 1/sqrt(10);
        [I,Q] = meshgrid([-3, -1, 1, 3], [-3, -1, 1, 3]);
        IQ_map = R*(I + 1j*Q);
        gray_idx = [0, 1, 3, 2, 4, 5, 7, 6, 12, 13, 15, 14, 8, 9, 11, 10];
        IQcodes(gray_idx+1) = IQ_map(:);
    elseif(isequal(modtype, '64QAM'))
        Nstates = 64; Nbits = 6; R = 1/sqrt(42);
        [I,Q] = meshgrid([-7,-5,-3,-1,1,3,5,7], [-7,-5,-3,-1,1,3,5,7]);
        IQ_map = R*(I + 1j*Q);
        IQcodes = IQ_map(:); % Non-gray
    elseif(isequal(modtype, '256QAM'))
        % --- Added for Exercise 20.15 ---
        Nstates = 256; Nbits = 8; R = 1/sqrt(170);
        % Levels are [-15, -13, ..., -1, 1, ..., 13, 15]
        [I,Q] = meshgrid(-15:2:15, -15:2:15);
        IQ_map = R*(I + 1j*Q);
        IQcodes = IQ_map(:); % Non-gray
    else
        error('Unknown modulation type: %s', modtype);
    end
end

function IQk = numbers2IQ(numbers, modtype, IQstates)
    if(isequal(modtype, 'DQPSK'))
        IQk(1) = 1; 
        for k = 1:length(numbers)
            IQk(k+1) = IQk(k) * IQstates(numbers(k)+1);
        end
    else
        IQk = IQstates(numbers + 1);
    end
    IQk = IQk(:);
end

function numbers = IQ2numbers(IQ, modtype, IQcodes)
    % This is the generic nearest-neighbor decoder.
    % It's robust and works for all modulations.
    N = length(IQ);
    numbers = zeros(N, 1);
    if size(IQcodes, 1) < size(IQcodes, 2), IQcodes = IQcodes.'; end

    % --- Specific (faster) decoders from book ---
    if(isequal(modtype, '16QAM'))
        R = 1/sqrt(10);
        levels = R * [-2, 0, 2]; % Decision boundaries
        % This is the corrected map from our debugging
        gray_map = [ 0  4 12  8;
                     1  5 13  9;
                     3  7 15 11;
                     2  6 14 10];
        for ns = 1:N
            I = real(IQ(ns)); Q = imag(IQ(ns));
            i_idx = sum(I > levels) + 1; % 1, 2, 3, or 4
            q_idx = sum(Q > levels) + 1; % 1, 2, 3, or 4
            numbers(ns) = gray_map(q_idx, i_idx);
        end
    elseif(isequal(modtype, '64QAM'))
        R = 1/sqrt(42);
        levels = R * [-6, -4, -2, 0, 2, 4, 6]; % 7 boundaries
        for ns = 1:N
            I = real(IQ(ns)); Q = imag(IQ(ns));
            i_idx = sum(I > levels); % 0-7
            q_idx = sum(Q > levels); % 0-7
            numbers(ns) = i_idx * 8 + q_idx; % Matches non-Gray map
        end
    elseif(isequal(modtype, '256QAM'))
        % --- Added for Exercise 20.15 ---
        R = 1/sqrt(170);
        levels = R * (-14:2:14); % 15 boundaries
        for ns = 1:N
            I = real(IQ(ns)); Q = imag(IQ(ns));
            i_idx = sum(I > levels); % 0-15
            q_idx = sum(Q > levels); % 0-15
            numbers(ns) = i_idx * 16 + q_idx; % Matches non-Gray map
        end
    else
        % Generic nearest-neighbor for other types (PAM, PSK, etc)
        for ns = 1:N
            [~, min_idx] = min(abs(IQ(ns) - IQcodes));
            numbers(ns) = min_idx - 1;
        end
    end
end

function plot_eye(signal, K, plot_title)
    % Helper function to plot an eye diagram
    % K = samples per symbol
    len = floor(length(signal) / (2*K)) * (2*K);
    if len == 0
        fprintf('Warning: Signal too short for eye diagram.\n');
        return;
    end
    sig_matrix = reshape(signal(1:len), 2*K, []);
    
    % Plot 2 symbol durations
    t = (0:2*K-1) / K; 
    plot(t, sig_matrix, 'b');
    title(plot_title);
    xlabel('Symbols');
    ylabel('Amplitude');
    grid on;
end