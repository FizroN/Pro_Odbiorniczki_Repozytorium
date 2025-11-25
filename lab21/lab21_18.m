% lab21_ex_equ_mse_nlms_robust_with_RLS.m
% Robust Channel Estimation and Equalization
% Methods: ZF, Regularized LS, LMS, NLMS, and **RLS**

% pkg load communications % Remove if using MATLAB
clear all; close all; clc;

% --- Parameters ---
N = 3000;    % Number of symbols
L = 1000;    % Training length for block equalizers
M = 10;      % Equalizer filter length                   
mi_nlms = 0.1;   % NLMS step size
mi_lms  = 0.005; % LMS step size 
delta   = 5.0;   % Damping factor for Regularized LS

% --- 1. Transmitter ---
x0 = (2*round(rand(1,N))-1) + j*(2*round(rand(1,N))-1); 
figure(1); subplot(2,2,1); plot(x0,'b*'); title('Input Constellation'); axis square;

% --- 2. Channel ---
h = [0.5; 0.2+0.2j; 0.1-0.05j; zeros(2,1); 0.05]; 
y_clean = conv(x0, h);

% --- 3. Receiver ---
y = awgn(y_clean, 20, 'measured'); 
subplot(2,2,2); plot(y,'r*'); title('Received (distorted + noisy)'); axis square;

% Prepare Data Matrices for Block Methods
V = toeplitz(y(M:L), y(M:-1:1));
u = x0(M:L).';

% --- METHOD A: Least Squares ---
g_zf = inv(V'*V) * V' * u; 
xest_zf = conv(y, g_zf.');

% --- METHOD B: Regularized LS ---
g_reg = inv(V'*V + delta * eye(M)) * V' * u;
xest_reg = conv(y, g_reg.');

subplot(2,2,3); plot(xest_zf(M:N), 'g.'); title('Standard LS (ZF)'); axis square;
subplot(2,2,4); plot(xest_reg(M:N), 'k.'); title('Regularized LS (Robust)'); axis square;

err_zf  = sum(abs(round(x0(M:N)) - round(xest_zf(M:N))) > 0.1);
err_reg = sum(abs(round(x0(M:N)) - round(xest_reg(M:N))) > 0.1);
fprintf('Bit Errors (Standard LS): %d\n', err_zf);
fprintf('Bit Errors (Regularized): %d\n', err_reg);
pause(0.1); 

% ==================================================================
% METHOD C: Adaptive Filters (LMS, NLMS, and RLS)
% ==================================================================
figure(2);

% Initialize Buffers
by = zeros(M,1);
xest_lms = zeros(1,N);
xest_nlms = zeros(1,N);
xest_rls  = zeros(1,N); % Output buffer for RLS

% Initialize Weights
g_lms = zeros(M,1);
g_nlms = zeros(M,1);
g_rls  = zeros(M,1);    % RLS Weights

% --- RLS Specific Parameters ---
P = 100 * eye(M);       % Inverse Correlation Matrix (Start large)
lambda = 0.99;          % Forgetting factor (0.98 to 1.0)

% Error history for convergence plot
mse_curve_nlms = zeros(1,N);
mse_curve_rls  = zeros(1,N);

for n = 1 : N
   by = [y(n); by(1:M-1)];  % Shift buffer
   
   % --- 1. Standard LMS ---
   y_out_lms = sum(by .* conj(g_lms));
   e_lms = x0(n) - y_out_lms;
   g_lms = g_lms + mi_lms * conj(e_lms) * by; 
   xest_lms(n) = y_out_lms;
   
   % --- 2. NLMS (Normalized) ---
   y_out_nlms = sum(by .* conj(g_nlms));
   e_nlms = x0(n) - y_out_nlms;
   energy = (by' * by) + 1e-5; 
   g_nlms = g_nlms + mi_nlms * conj(e_nlms) * by / energy; 
   xest_nlms(n) = y_out_nlms;
   
   mse_curve_nlms(n) = abs(e_nlms)^2; % Save error for plot

   % --- 3. RLS (Recursive Least Squares) [NEW PART] ---
   % A. Filter Operation
   y_out_rls = sum(by .* conj(g_rls));
   e_rls = x0(n) - y_out_rls;
   xest_rls(n) = y_out_rls;
   
   mse_curve_rls(n) = abs(e_rls)^2; % Save error for plot
   
   % B. Calculate Gain Vector (k)
   % This "whitens" the input using the P matrix history
   Pi = P * by;
   k = Pi / (lambda + by' * Pi);
   
   % C. Update Weights
   g_rls = g_rls + k * conj(e_rls);
   
   % D. Update Inverse Correlation Matrix P
   P = (P - k * by' * P) / lambda;
   
end

% Plotting Adaptive Results
subplot(3,1,1); plot(xest_lms(N/2:end), 'm.'); 
title('LMS Output'); axis([-2 2 -2 2]); axis square;

subplot(3,1,2); plot(xest_nlms(N/2:end), 'b.'); 
title('NLMS Output'); axis([-2 2 -2 2]); axis square;

subplot(3,1,3); plot(xest_rls(N/2:end), 'r.'); 
title('RLS Output (Note tighter clusters)'); axis([-2 2 -2 2]); axis square;

% Calculate final errors
err_lms = sum(round(x0(N/2:end)) ~= round(xest_lms(N/2:end)));
err_nlms = sum(round(x0(N/2:end)) ~= round(xest_nlms(N/2:end)));
err_rls = sum(round(x0(N/2:end)) ~= round(xest_rls(N/2:end)));

fprintf('Errors LMS:  %d\n', err_lms);
fprintf('Errors NLMS: %d\n', err_nlms);
fprintf('Errors RLS:  %d  <-- Should be lowest/similar to NLMS but faster\n', err_rls);

% ==================================================================
% NEW: Convergence Speed Comparison Plot
% ==================================================================
figure(3);
% Smooth the error curve to make it readable
window = 20; 
smooth_nlms = filter(ones(1,window)/window, 1, mse_curve_nlms);
smooth_rls  = filter(ones(1,window)/window, 1, mse_curve_rls);

semilogy(smooth_nlms, 'b', 'LineWidth', 1); hold on;
semilogy(smooth_rls, 'r', 'LineWidth', 2);
grid on;
legend('NLMS (Gradient Descent)', 'RLS (Matrix Inverse)');
title('Convergence Speed Comparison');
xlabel('Symbol Iteration'); ylabel('Mean Squared Error (dB)');
axis([0 500 10^-3 10^1]); % Zoom in on the first 500 symbols

fprintf('\nObserve Figure 3: RLS (Red) drops to the noise floor much faster than NLMS (Blue).\n');