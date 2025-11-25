% lab21_ex_equ_mse_nlms_robust.m
% Robust Channel Estimation and Equalization
% Methods: ZF (via LS), Regularized LS (Damped), LMS, NLMS

%pkg load communications % Remove if using MATLAB
clear all; close all;

% --- Parameters ---
N = 3000;    % Number of symbols
L = 1000;    % Training length for block equalizers
M = 10;      % Equalizer filter length                   
mi_nlms = 0.1;   % NLMS step size
mi_lms  = 0.005; % LMS step size (must be smaller than NLMS)
delta   = 5.0;   % Damping factor for Regularized LS (Robustness)

% --- 1. Transmitter ---
% IQ 4-QAM states
x0 = (2*round(rand(1,N))-1) + j*(2*round(rand(1,N))-1); 
figure(1); subplot(2,2,1); plot(x0,'b*'); title('Input Constellation'); axis square;

% --- 2. Channel (Multipath) ---
h = [0.5; 0.2+0.2j; 0.1-0.05j; zeros(2,1); 0.05]; % Fixed complex channel
y_clean = conv(x0, h);

% --- 3. Receiver (Add Noise) ---
% SNR is lower here (20dB) to highlight the benefit of robust algorithms
y = awgn(y_clean, 20, 'measured'); 
subplot(2,2,2); plot(y,'r*'); title('Received (distorted + noisy)'); axis square;

% Prepare Data Matrices for Block Methods
V = toeplitz(y(M:L), y(M:-1:1));  % Convolution matrix (received)
u = x0(M:L).';                    % Target vector (transmitted)

% ==================================================================
% METHOD A: Least Squares (Zero Forcing approximation)
% ==================================================================
% This is what V\u does. It tries to force error to zero, 
% but amplifies noise if matrix is ill-conditioned.
g_zf = inv(V'*V) * V' * u; 
xest_zf = conv(y, g_zf.');

% ==================================================================
% METHOD B: Regularized Least Squares (Robust / "Damped")
% ==================================================================
% We add 'delta * Identity Matrix' to the inversion.
% This prevents the filter from creating huge gains for noise.
% Formula: g = inv(V'*V + delta*I) * V'*u
g_reg = inv(V'*V + delta * eye(M)) * V' * u;
xest_reg = conv(y, g_reg.');

subplot(2,2,3); plot(xest_zf(M:N), 'g.'); title('Standard LS (ZF)'); axis square;
subplot(2,2,4); plot(xest_reg(M:N), 'k.'); title('Regularized LS (Robust)'); axis square;

% Compare Errors
err_zf  = sum(abs(round(x0(M:N)) - round(xest_zf(M:N))) > 0.1);
err_reg = sum(abs(round(x0(M:N)) - round(xest_reg(M:N))) > 0.1);
fprintf('Bit Errors (Standard LS): %d\n', err_zf);
fprintf('Bit Errors (Regularized): %d\n', err_reg);
fprintf('--> Regularization usually reduces outliers.\n\n');

pause; 

% ==================================================================
% METHOD C: Adaptive Filters (LMS vs NLMS)
% ==================================================================
figure(2);

% Initialize
by = zeros(M,1);
g_lms = zeros(M,1);
g_nlms = zeros(M,1);
xest_lms = zeros(1,N);
xest_nlms = zeros(1,N);

for n = 1 : N
   by = [y(n); by(1:M-1)];  % Shift buffer
   
   % --- 1. Standard LMS ---
   % Simple gradient descent. 
   % Update: w = w + mu * error * input
   y_out_lms = sum(by .* conj(g_lms));
   e_lms = x0(n) - y_out_lms;
   g_lms = g_lms + mi_lms * conj(e_lms) * by; 
   xest_lms(n) = y_out_lms;
   
   % --- 2. NLMS (Normalized) ---
   % Robust to input power scaling.
   % Update: w = w + mu * error * input / (input_energy)
   y_out_nlms = sum(by .* conj(g_nlms));
   e_nlms = x0(n) - y_out_nlms;
   energy = (by' * by) + 1e-5; % Small constant to avoid div by zero
   g_nlms = g_nlms + mi_nlms * conj(e_nlms) * by / energy; 
   xest_nlms(n) = y_out_nlms;
end

subplot(2,1,1); plot(xest_lms(N/2:end), 'm.'); 
title('LMS Equalizer Output (Converged)'); axis([-2 2 -2 2]); axis square;
subplot(2,1,2); plot(xest_nlms(N/2:end), 'b.'); 
title('NLMS Equalizer Output (Converged)'); axis([-2 2 -2 2]); axis square;

% Calculate final errors
err_lms = sum(round(x0(N/2:end)) ~= round(xest_lms(N/2:end)));
err_nlms = sum(round(x0(N/2:end)) ~= round(xest_nlms(N/2:end)));

fprintf('Errors LMS: %d\n', err_lms);
fprintf('Errors NLMS: %d\n', err_nlms);