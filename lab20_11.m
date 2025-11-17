% lab20_ex_adc.m
% This script models the effect of an ADC sampling rate error 
% and then demonstrates its perfect correction.
% This is the analysis for Exercise 20.11.

clear all; close all;

% --- 1. Parameters ---
N = 1000; n = 0:N-1; % Time vector, 1000 samples
fs = 1e+6;           % Ideal sampling frequency (1 MHz)
fc = 1e+5;           % Ideal carrier frequency (100 kHz)
ADCppm = 100;        % The ADC error: 100 parts-per-million

% --- 2. The Model ---
% 'A' is the epsilon (e) in Eq (20.26)
% A = 100 * 1e-6 = 0.0001
A = ADCppm * 1e-6;   

% 'x1' is the PERFECT, ideal signal we *wanted* to receive.
% The frequency is exactly fc/fs.
x1 = exp(j*2*pi*fc/fs*n);

% 'x2' is the BAD signal we *actually* received.
% The ADC error 'A' means the *actual* sampling rate is (1+A)*fs.
% This is y1(n) from Eq (20.23).
x2 = exp(j*2*pi*fc/((1+A)*fs)*n);

% --- 3. The Correction ---
% From Eq (20.24) and (20.25), we can correct 'x2' by multiplying
% it by a complex sinusoid with the *opposite* frequency offset.
% The frequency offset 'x' in the text is (A/(1+A))*(fc/fs).
% (Note: text has a minus sign, but we apply the opposite for correction).
correction_factor = exp(j*2*pi * A/(1+A)*fc/fs * n);
x1c = x2 .* correction_factor; % This is the corrected signal

% --- 4. The Proof ---
% If the correction is perfect, the difference between the
% ideal signal (x1) and the corrected signal (x1c) should be zero.
error = max(abs(x1-x1c));

fprintf('Ideal sampling freq (fs): %.0f Hz\n', fs);
fprintf('Actual sampling freq ((1+A)*fs): %.0f Hz\n', (1+A)*fs);
fprintf('ADC Error (A): %f\n', A);
fprintf('Max error between ideal and corrected signal: %e\n', error);

% The error will be a very small number (e.g., 1.4e-14),
% which is just machine rounding error. This proves the correction works.