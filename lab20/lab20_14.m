% lab20_ex_image_tx.m
% This script is for Exercise 20.14: Image transmission.
% It modifies the full transmission chain to send an image file
% (the built-in MATLAB 'peppers.png' image) instead of text or audio.

clear all; close all;

% =====================================================================
% 1. Load Image Payload
% =====================================================================
fprintf('Loading image payload...\n');
% Load a built-in MATLAB image
try
    img_in = imread('peppers.png');
catch
    fprintf('Could not find peppers.png, using random image.\n');
    img_in = randi(255, [256, 256, 3], 'uint8');
end

% Get image dimensions to reconstruct it later
[img_height, img_width, img_channels] = size(img_in);
fprintf('Loaded image: %d x %d x %d\n', img_height, img_width, img_channels);

% Serialize the 3D matrix into a 1D stream of bytes
% This is the data we will transmit
payload_bytes_in = img_in(:); % Flattens the matrix
payload_len_bytes = length(payload_bytes_in);

fprintf('Serialized image into %d bytes.\n', payload_len_bytes);

% =====================================================================
% 2. Simulation Parameters
% =====================================================================
% --- Modulation Parameters ---
modtype = '16QAM'; % 4QAM, 16QAM, 64QAM. Higher QAM = more bits/symbol

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
[IQcodes, Nstates, Nbits, R] = IQdef(modtype);

% --- Convert payload (bytes) into modulation symbols (numbers) ---
[numbers, ~] = bytes2numbers(payload_bytes_in, Nbits);
fprintf('Packed %d bytes into %d symbols using %d-bits/symbol (%s).\n', ...
        payload_len_bytes, length(numbers), Nbits, modtype);

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

% Decode IQ symbols back to numbers
received_numbers = IQ2numbers(IQkk_received, modtype, IQcodes);

% Convert numbers back to audio bytes
payload_bytes_out = numbers2bytes(received_numbers, Nbits, payload_len_bytes);

fprintf('Received and decoded %d bytes.\n', length(payload_bytes_out));

% =====================================================================
% 6. Verification (View Images)
% =====================================================================

% --- Check if we received the right number of bytes ---
if length(payload_bytes_out) < payload_len_bytes
    fprintf('Error: Received fewer bytes than expected!\n');
    % Pad with zeros if reception was incomplete
    payload_bytes_out(payload_len_bytes) = 0;
end
% Trim any extra padding bytes
payload_bytes_out_trimmed = payload_bytes_out(1:payload_len_bytes);


% --- Deserialize: Reshape the 1D byte stream back into the 3D image ---
img_out = reshape(payload_bytes_out_trimmed, img_height, img_width, img_channels);

% --- Plot results ---
figure;
subplot(1,2,1);
imshow(img_in);
title('Original Image');

subplot(1,2,2);
imshow(img_out);
title(sprintf('Received Image (SNR: %d dB, Mod: %s)', SNR_dB, modtype));

fprintf('Displaying original and received images.\n');
fprintf('Done.\n');


% =====================================================================
% HELPER FUNCTIONS (Same as lab20_ex_audio_tx.m)
% =====================================================================

function [numbers, bitsnum] = bytes2numbers(bytes_in, Nbits)
    % --- FIX START ---
    % The previous logic was interlacing the bits (bit-plane-wise).
    % This is the correct logic to flatten the bytes sequentially
    % (byte-wise).
    
    % Convert array of bytes to a char matrix, e.g.:
    % '01001000' (byte 1)
    % '01100101' (byte 2)
    % ...and transpose it to [8, N_bytes]
    bits_char_matrix = dec2bin(bytes_in, 8)';
    
    % Reshape into a single long stream of '0' and '1' chars
    % This reads column 1 (byte 1), then column 2 (byte 2), etc.
    bits_stream = reshape(bits_char_matrix, [1, numel(bits_char_matrix)]);
    % --- FIX END ---
    
    N = length(bits_stream);
    Nadd = Nbits - rem(N, Nbits);
    if Nadd == Nbits, Nadd = 0; end
    for k = 1:Nadd, bits_stream = [bits_stream, '0']; end
    bitsnum_matrix = reshape(bits_stream, [Nbits, length(bits_stream)/Nbits])';
    numbers = bin2dec(bitsnum_matrix);
    bitsnum = bitsnum_matrix; % for debugging
end

function bytes_out = numbers2bytes(numbers, Nbits, expected_byte_len)
    text_bits = dec2bin(numbers, Nbits);
    [rows, cols] = size(text_bits);

    % --- FIX START ---
    % Original code used ' and reshape, which shuffled bits
    % Transpose, then flatten with (:) to get the correct 1D bit stream
    bits_stream_matrix_t = text_bits'; % [Nbits, N_syms]
    bits_stream = bits_stream_matrix_t(:)'; % 1D stream
    % --- FIX END ---
    
    N_bytes = min(expected_byte_len, floor(length(bits_stream) / 8));
    N_bits = N_bytes * 8;
    if N_bits == 0
        bytes_out = uint8([]);
        return;
    end
    bits_stream = bits_stream(1:N_bits);
    
    % --- FIX START ---
    % Reshape into an [8, N_bytes] matrix, then transpose
    bits_stream_matrix = reshape(bits_stream, [8, N_bytes]); % [8, N_bytes]
    text_bytes_matrix = bits_stream_matrix'; % [N_bytes, 8]
    % --- FIX END ---
    
    bytes_out = uint8(bin2dec(text_bytes_matrix));
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
    N = length(IQ);
    numbers = zeros(N, 1);
    if size(IQcodes, 1) < size(IQcodes, 2), IQcodes = IQcodes.'; end

    % Use specific decoders where available, otherwise generic
    if(isequal(modtype, '16QAM'))
        R = 1/sqrt(10);
        levels = R * [-2, 0, 2]; % Decision boundaries
        
        % --- FIX START ---
        % The original gray_map was the transpose of what it should be,
        % causing a mismatch with the encoder's mapping.
        % This is the correct map, based on the encoder's gray_idx
        gray_map = [ 0  4 12  8;
                     1  5 13  9;
                     3  7 15 11;
                     2  6 14 10];
        % --- FIX END ---
                     
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
    else
        % Generic nearest-neighbor for other types (PAM, PSK, etc)
        for ns = 1:N
            [~, min_idx] = min(abs(IQ(ns) - IQcodes));
            numbers(ns) = min_idx - 1;
        end
    end
end