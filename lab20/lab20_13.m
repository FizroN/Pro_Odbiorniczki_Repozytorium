% lab20_ex_audio_tx.m
% This script is for Exercise 20.13: Speech/audio transmission.
% It modifies the full transmission chain to send an audio file
% (the built-in MATLAB 'chirp' sound) instead of text.

clear all; close all;

% =====================================================================
% 1. Load Audio Payload
% =====================================================================
fprintf('Loading audio payload...\n');
% Load a built-in MATLAB sound
load chirp; 
payload_Fs = Fs;
% Take the first 8192 samples (approx 1 sec)
payload_len = 8192; 
payload_in = y(1:payload_len);
% Quantize: Scale double [-1, 1] to uint8 [0, 255]
payload_bytes_in = uint8((payload_in + 1) / 2 * 255);

fprintf('Loaded %d audio samples (bytes).\n', length(payload_bytes_in));

% =====================================================================
% 2. Simulation Parameters (from lab20_ex_pulse_shaping)
% =====================================================================
% --- Modulation Parameters ---
modtype = '4QAM'; % 2PAM, 4PAM, 8PAM, BPSK, QPSK, 8PSK, 4QAM, 16QAM

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
SNR_dB = 160; chan_gain = 1; chan_phase = 0; carDf = 0; carDph = 0; ADCppm = 0;

% Disturbed values (uncomment to test "rainy days"):
%SNR_dB = 20; chan_gain = 0.25; chan_phase = pi/7; carDf = 0.001; carDph = pi/11; ADCppm = 100;
fprintf('--- Using Disturbances ---\n');
fprintf('SNR: %d dB\n', SNR_dB);

% =====================================================================
% 3. Transmitter Side
% =====================================================================
fprintf('Transmitting...\n');
% Get modulation definitions
[IQcodes, Nstates, Nbits, R] = IQdef(modtype);

% --- Convert payload (bytes) into modulation symbols (numbers) ---
% This is the modified 'text2numbers'
[numbers, bits_packed] = bytes2numbers(payload_bytes_in, Nbits);
fprintf('Packed %d bytes into %d symbols using %d-bits/symbol (%s).\n', ...
        length(payload_bytes_in), length(numbers), Nbits, modtype);

% Generate IQ symbols
IQk = numbers2IQ(numbers, modtype, IQcodes);

% --- Pulse Shaping (from lab20_ex_pulse_shaping) ---
fcut = fs/(2*K); % PSF cut-off frequency
Npsf = Ns*K+1; Mpsf = (Npsf-1)/2; % PSF filter length and its half

% Add prefix/postfix to handle filter ramp-up/down
N_prefix = Ns;
N_postfix = Ns + ceil( (N_prefix*K + length(IQk)*K + N_prefix*K) / (N_prefix*K) );
numbers_dummy = floor(Nstates*(rand(2*N_prefix,1)-10*eps));
dummy = numbers2IQ(numbers_dummy, modtype, IQcodes);
IQdum = [dummy(1:N_prefix); IQk; dummy(N_prefix+1:end)]; % dummy prefix and postfix

IQ0 = zeros(1, length(IQdum)*K); % appending prefix & postfix
IQ0(1:K:end) = IQdum; % zero insertion

% Generate Filter
if(isequal(psf_type, 'rect'))
    hpsf = [ones(1,K) zeros(1, (Ns-1)*K+1)] / K;
    hpsf(Npsf) = 0;
else
    rctype = 'sqrt';
    if isequal(psf_type, 'normal-rc'), rctype = 'normal'; end
    hpsf = firrcos(Npsf-1, fcut, r, fs, 'rolloff', rctype);
end

% Apply filter
IQn = conv(IQ0, hpsf); 
IQn = IQn(Mpsf+1 : end-Mpsf); % pulse shaping in TX

% =====================================================================
% 4. Channel Simulation (from lab20_ex_pulse_shaping)
% =====================================================================

if (do_updown) 
    fprintf('Simulating Frequency UP/DOWN conversion...\n');
    N = length(IQn); n = 0:N-1;
    
    % --- Frequency UP conversion (TX) ---
    y = real(IQn).*cos(2*pi*fcar/fs*n) - imag(IQn).*sin(2*pi*fcar/fs*n);

    % --- Frequency DOWN conversion (RX) ---
    df = 0; dphi = 0; 
    IQnn_raw = 2*y.*exp(-j*(2*pi*(fcar/fs + df/fs)*n + dphi));
    IQnn = IQnn_raw;
    
    if(do_disturb)
        fprintf('Applying disturbances...\n');
        % 1. Channel
        IQnn = IQnn * chan_gain * exp(j*chan_phase);
        % 2. Noise
        s_power = mean(abs(IQnn).^2);
        if s_power == 0, s_power = 1; end 
        noise_power = s_power / (10^(SNR_dB/10));
        noise = sqrt(noise_power/2) * (randn(size(IQnn)) + 1j*randn(size(IQnn)));
        IQnn = IQnn + noise;
        % 3. Receiver Errors
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
% Apply matched filter (the same PSF)
IQnn = conv(IQnn, hpsf); 
IQnn = IQnn(Mpsf+1 : end-Mpsf); % pulse shaping in RX

% Synchronization and Sampling (ideal)
IQnn_synced = IQnn(N_prefix*K+1 : end-N_prefix*K); % Remove ramp-up/down
IQkk_received = IQnn_synced(1:K:end); % Sample at symbol rate

% --- Decode IQ symbols back to numbers ---
received_numbers = IQ2numbers(IQkk_received, modtype, IQcodes);

% --- Convert numbers back to audio bytes ---
% This is the modified 'numbers2text'
payload_bytes_out = numbers2bytes(received_numbers, Nbits, length(payload_bytes_in));

fprintf('Received and decoded %d bytes.\n', length(payload_bytes_out));

% =====================================================================
% 6. Verification (Listen and Plot)
% =====================================================================

% --- Convert uint8 [0, 255] back to double [-1, 1] ---
payload_out = (double(payload_bytes_out) / 255) * 2 - 1;

% --- Plot results ---
figure;
subplot(2,1,1);
plot(payload_in, 'b');
title(sprintf('Original Audio (First %d samples)', payload_len));
xlabel('Sample'); ylabel('Amplitude');
ylim([-1.2 1.2]); grid on;

subplot(2,1,2);
plot(payload_out, 'r');
title(sprintf('Received Audio (SNR: %d dB, Mod: %s)', SNR_dB, modtype));
xlabel('Sample'); ylabel('Amplitude');
ylim([-1.2 1.2]); grid on;

% --- Play sounds ---
fprintf('Playing original sound...\n');
sound(payload_in, payload_Fs);
pause(length(payload_in)/payload_Fs + 0.5); % Wait for sound to finish

fprintf('Playing received sound...\n');
sound(payload_out, payload_Fs);


% =====================================================================
% HELPER FUNCTIONS
% =====================================================================

function [numbers, bitsnum] = bytes2numbers(bytes_in, Nbits)
    % This function converts an array of uint8 bytes into
    % an array of 'numbers' (symbols), where each symbol
    % represents Nbits.
    
    % Convert array of bytes [72, 101, 108] to a matrix of bit-strings
    % '01001000'
    % '01100101'
    % '01101100'
    bits_char_matrix = dec2bin(bytes_in, 8)';
    
    % Reshape into a single long stream of '0' and '1' chars
    % '010010000110010101101100...'
    bits_stream = reshape(bits_char_matrix, [1, numel(bits_char_matrix)]);
    
    % Pad with '0's at the end to make length a multiple of Nbits
    N = length(bits_stream);
    Nadd = Nbits - rem(N, Nbits);
    if Nadd == Nbits, Nadd = 0; end % Handle case where N is multiple of Nbits
    for k = 1:Nadd, bits_stream = [bits_stream, '0']; end
    
    % Reshape into a matrix where each column is Nbits
    % e.g., for Nbits=4:
    % '0100'
    % '1000'
    % '0110' ...
    bitsnum_matrix = reshape(bits_stream, [Nbits, length(bits_stream)/Nbits])';
    
    % Convert binary strings to decimal numbers
    numbers = bin2dec(bitsnum_matrix);
    bitsnum = bitsnum_matrix; % for debugging
end

function bytes_out = numbers2bytes(numbers, Nbits, expected_byte_len)
    % This is the reverse of bytes2numbers.
    % It converts an array of symbols ('numbers') back into
    % an array of uint8 bytes.
    
    % Convert decimal numbers [4, 8, 6] to binary strings
    % '0100'
    % '1000'
    % '0110'
    text_bits = dec2bin(numbers, Nbits);
    
    % Reshape into a single long stream of '0' and '1' chars
    [rows, cols] = size(text_bits);
    bits_stream = reshape(text_bits', [1, rows*cols]);
    
    % We need to cut it to the nearest 8-bit boundary,
    % and not exceed the original number of *bits*.
    N_bytes = min(expected_byte_len, floor(length(bits_stream) / 8));
    N_bits = N_bytes * 8;
    
    % Trim to the last full byte
    bits_stream = bits_stream(1:N_bits);
    
    % Reshape into 8-bit rows
    % '01001000'
    % '01100101'
    % '01101100'
    text_bytes_matrix = reshape(bits_stream, [8, N_bytes])';
    
    % Convert binary strings to decimal numbers and cast to uint8
    bytes_out = uint8(bin2dec(text_bytes_matrix));
end


% --- Other helper functions (IQdef, numbers2IQ, IQ2numbers) ---

function [IQcodes, Nstates, Nbits, R] = IQdef(modtype)
    % Defines the constellation points (IQ codes) for modulations
    Nstates = 0; Nbits = 0; R = 1; IQcodes = [];
    if(isequal(modtype, '2PAM'))
        Nstates = 2; Nbits = 1; IQcodes = [-1, 1];
    elseif(isequal(modtype, '4PAM'))
        Nstates = 4; Nbits = 2; IQcodes = [-3, -1, 3, 1]; % Gray code
    elseif(isequal(modtype, '8PAM'))
        Nstates = 8; Nbits = 3; IQcodes = [-7, -5, -3, -1, 7, 5, 3, 1]; % Gray
    elseif(isequal(modtype, 'BPSK'))
        Nstates = 2; Nbits = 1; IQcodes = [1, -1]; % 0, pi
    elseif(isequal(modtype, 'QPSK'))
        Nstates = 4; Nbits = 2; R = 1/sqrt(2); 
        IQcodes = [R+j*R, -R+j*R, -R-j*R, R-j*R]; % R*exp(j*pi/4*[1 3 5 7])
    elseif(isequal(modtype, 'DQPSK'))
        Nstates = 4; Nbits = 2; R = 1;
        IQcodes = [exp(j*pi/4), exp(j*3*pi/4), exp(j*7*pi/4), exp(j*5*pi/4)];
    elseif(isequal(modtype, '8PSK'))
        Nstates = 8; Nbits = 3; 
        idx = [0, 1, 3, 2, 7, 6, 4, 5]; % Gray code
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
    else
        error('Unknown modulation type: %s', modtype);
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

function numbers = IQ2numbers(IQ, modtype, IQcodes)
    % from [I,Q] values to carrier state numbers for many input IQ pairs
    % This function finds the *closest* constellation point.
    
    N = length(IQ);
    numbers = zeros(N, 1);
    
    % A generic (but slower) nearest-neighbor solution
    % This works for all modulation types defined in IQdef
    
    % Re-center IQcodes to ensure it's a column vector for distance calc
    if size(IQcodes, 1) < size(IQcodes, 2)
        IQcodes = IQcodes.';
    end

    for ns = 1:N
        % Calculate Euclidean distance from the received symbol
        % to all possible ideal constellation points
        distances = abs(IQ(ns) - IQcodes);
        
        % Find the index of the minimum distance
        [~, min_idx] = min(distances);
        
        % The number is the index (0-based)
        numbers(ns) = min_idx - 1;
    end
    
    % --- Specific (faster) decoders from book ---
    % Note: The generic solution above is more robust, but if we
    % want to use the exact (and faster) book logic, we can:
    
    if(isequal(modtype, '16QAM'))
        R = 1/sqrt(10);
        levels = R * [-2, 0, 2]; % Decision boundaries
        gray_map = [ 0  1  3  2;
                     4  5  7  6;
                    12 13 15 14;
                     8  9 11 10];
        for ns = 1:N
            I = real(IQ(ns));
            Q = imag(IQ(ns));
            
            i_idx = sum(I > levels) + 1; % 1, 2, 3, or 4
            q_idx = sum(Q > levels) + 1; % 1, 2, 3, or 4
            
            numbers(ns) = gray_map(q_idx, i_idx);
        end
    elseif(isequal(modtype, '64QAM'))
        R = 1/sqrt(42);
        levels = R * [-6, -4, -2, 0, 2, 4, 6]; % 7 boundaries
        for ns = 1:N
            I = real(IQ(ns));
            Q = imag(IQ(ns));
            
            i_idx = sum(I > levels); % 0-7
            q_idx = sum(Q > levels); % 0-7

            % This must match the (non-Gray) mapping in IQdef
            % which was IQ_map = R*(I + 1j*Q); IQcodes = IQ_map(:);
            % meshgrid(I,Q) makes I vary by column, Q by row
            % So IQ_map(:) is [Icol1, Icol2, ...]
            numbers(ns) = i_idx * 8 + q_idx;
        end
    end
end