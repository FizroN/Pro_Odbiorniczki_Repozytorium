% lab20_ex_fsk.m
% This script is for Exercise 20.16: Frequency Shift Keying (FSK)
% It demonstrates the core concept of FSK by changing the
% *rate of phase change* to encode bits.

clear all; close all;

% =====================================================================
% 1. Parameters
% =====================================================================
fprintf('Setting up FSK parameters...\n');

bits_in = [1, 0, 1, 1, 0, 1, 0, 0, 1]; % Bitstream to transmit
K = 50; % Samples per symbol (per bit)
fs = 1000; % Sampling frequency (Hz)

f_carrier = 100; % A central carrier frequency (Hz)
f_dev = 30;      % Frequency deviation (Hz)

% Define the two frequencies for '0' and '1'
f0 = f_carrier - f_dev; % Frequency for bit '0'
f1 = f_carrier + f_dev; % Frequency for bit '1'

fprintf('Bitstream: %s\n', num2str(bits_in));
fprintf('Samples per bit (K): %d\n', K);
fprintf('Frequency for 0: %.0f Hz\n', f0);
fprintf('Frequency for 1: %.0f Hz\n', f1);

% =====================================================================
% 2. FSK Modulator (Transmitter)
% =====================================================================
fprintf('Modulating...\n');

% Create the signal, one bit at a time
signal_out = [];
current_phase = 0; % We must maintain phase continuity

% Normalized angular frequencies
w0 = 2 * pi * f0 / fs;
w1 = 2 * pi * f1 / fs;

for bit = bits_in
    % Create time vector for this one symbol
    t = (0:K-1)';
    
    if bit == 0
        % Bit '0': Use frequency w0
        % Phase = w0*t + previous_phase
        symbol_phase = w0 * t + current_phase;
    else
        % Bit '1': Use frequency w1
        % Phase = w1*t + previous_phase
        symbol_phase = w1 * t + current_phase;
    end
    
    % Generate the cosine for this symbol
    symbol = cos(symbol_phase);
    
    % Add to our signal
    signal_out = [signal_out; symbol];
    
    % Store the *last* phase to ensure continuity
    current_phase = symbol_phase(end);
end

% --- Add noise (like in our other exercises) ---
SNR_dB = 10;
s_power = mean(abs(signal_out).^2);
noise_power = s_power / (10^(SNR_dB/10));
noise = sqrt(noise_power) * randn(size(signal_out));
signal_in = signal_out + noise;


% =====================================================================
% 3. FSK Demodulator (Receiver) - Simple (non-coherent)
% =====================================================================
% This is a simple (and not very robust) way to decode.
% We will use two bandpass filters, one centered at f0 and one at f1.
% Whichever filter has more energy, that's our bit.
fprintf('Demodulating...\n');

% Design filters
bp_0 = fir1(100, [f0-10, f0+10]/(fs/2), 'bandpass');
bp_1 = fir1(100, [f1-10, f1+10]/(fs/2), 'bandpass');

% Reshape received signal into a matrix (one symbol per column)
% We have to ignore the filter delay, so this is imperfect.
num_bits = length(bits_in);
signal_matrix = reshape(signal_in(1:num_bits*K), K, num_bits);

bits_out = [];
for i = 1:num_bits
    symbol = signal_matrix(:, i);
    
    % Filter the symbol with both filters
    energy_0 = sum(abs(filter(bp_0, 1, symbol)).^2);
    energy_1 = sum(abs(filter(bp_1, 1, symbol)).^2);
    
    % Make a decision
    if energy_1 > energy_0
        bits_out(i) = 1;
    else
        bits_out(i) = 0;
    end
end


% =====================================================================
% 4. Verification
% =====================================================================

% --- Plot modulated signal ---
figure;
subplot(2,1,1);
t_total = (0:length(signal_out)-1) / fs;
plot(t_total, signal_out);
hold on;
% Add bit boundaries
for i = 1:length(bits_in)
    plot([i*K/fs, i*K/fs], [-1.5 1.5], 'r:');
end
title(sprintf('Transmitted FSK Signal (f0=%.0f, f1=%.0f)', f0, f1));
xlabel('Time (s)'); ylabel('Amplitude');
ylim([-1.5 1.5]);
% Show the bits on the plot
text(0.01, 1.2, 'Bits:','FontSize',10);
for i = 1:length(bits_in)
    text((i*K - K/2)/fs, 1.2, num2str(bits_in(i)), 'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold');
end

% --- Plot received (noisy) signal ---
subplot(2,1,2);
plot(t_total, signal_in);
title(sprintf('Received FSK Signal (SNR: %d dB)', SNR_dB));
xlabel('Time (s)'); ylabel('Amplitude');
ylim([-3 3]); % Zoom out to see noise


% --- Print results ---
fprintf('\n--- TRANSMISSION COMPLETE ---\n');
fprintf('Original Bits: %s\n', num2str(bits_in));
fprintf('Received Bits: %s\n', num2str(bits_out));
num_errors = sum(bits_in ~= bits_out);
fprintf('Bit Errors: %d / %d\n', num_errors, num_bits);

fprintf('All exercises complete. Good work!\n');