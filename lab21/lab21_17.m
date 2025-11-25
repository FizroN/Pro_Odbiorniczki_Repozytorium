% lab21_complete_equalization.m
% Comparison of Channel Equalization Algorithms:
% 1. Block LS (Zero Forcing)
% 2. Block Regularized LS (Robust)
% 3. LMS (Least Mean Squares)
% 4. NLMS (Normalized LMS)
% 5. RLS (Recursive Least Squares) - NEW

pkg load communications % Only needed for Octave, comment out for MATLAB
clear all; close all; clc;

% --- Parameters ---
N = 4000;          % Number of symbols
M = 10;            % Filter length (Taps)
L_train = 500;     % Number of symbols used for Block Training

% Adaptive Parameters
mu_lms  = 0.005;   % LMS Step size (keep small for stability)
mu_nlms = 0.1;     % NLMS Step size (can be larger than LMS)
lambda_rls = 0.99; % RLS Forgetting factor (0.98 to 0.999 usually)

% --- 1. Signal Generation (QPSK/4-QAM) ---
% Create random +1/-1 symbols for I and Q
x0 = (2*round(rand(N,1))-1) + 1j*(2*round(rand(N,1))-1);

% --- 2. Channel Definition (Multipath) ---
% A static channel with some echoes
h = [0.5; 0.2+0.2j; 0.1-0.05j; zeros(2,1); 0.05]; 
y_clean = filter(h, 1, x0); % Apply channel

% --- 3. Add Noise ---
SNR_dB = 20;
y = awgn(y_clean, SNR_dB, 'measured');

% Visualization of input vs received
figure(1);
subplot(2,1,1); plot(x0(1:100), 'bo'); title('Transmitted Symbols'); axis square;
subplot(2,1,2); plot(y(1:1000), 'r.'); title('Received (Multipath + Noise)'); axis square;
pause(0.1);

% =================================================================
% PART A: Block-Based Methods (Offline Processing)
% =================================================================
disp('Running Block Methods...');

% Create Toeplitz Matrix using training data L_train
% We use the received data 'y' to try to recover transmitted 'x0'
recv_vector = y(M:L_train);
% Construct convolution matrix
V = toeplitz(recv_vector, y(M:-1:1)); 
desired = x0(M:L_train);

% 1. Standard Least Squares (Zero Forcing)
% g = (V'V)^-1 V'u
g_zf = (V'*V) \ (V' * desired);

% 2. Regularized Least Squares (Robust)
% g = (V'V + delta*I)^-1 V'u
delta = 5.0; % Damping factor
g_reg = inv(V'*V + delta*eye(M)) * V' * desired;

% Apply to whole signal
xest_zf  = filter(g_zf, 1, y);
xest_reg = filter(g_reg, 1, y);

% =================================================================
% PART B: Adaptive Methods (Online/Loop Processing)
% =================================================================
disp('Running Adaptive Methods (LMS, NLMS, RLS)...');

% Initialization
by = zeros(M,1);          % Input buffer (regressor)

% LMS Init
w_lms = zeros(M,1);
xest_lms = zeros(N,1);
err_lms_curve = zeros(N,1);

% NLMS Init
w_nlms = zeros(M,1);
xest_nlms = zeros(N,1);
err_nlms_curve = zeros(N,1);

% RLS Init
w_rls = zeros(M,1);
xest_rls = zeros(N,1);
P = 100 * eye(M);         % Inverse correlation matrix (start large)
err_rls_curve = zeros(N,1);

for n = 1:N
    % 1. Update Buffer (Slide window of received samples)
    % If n < M, pad with zeros, otherwise take recent M samples
    if n < M
        by = [y(n:-1:1); zeros(M-n,1)];
    else
        by = y(n:-1:n-M+1);
    end
    
    desired_sym = x0(n);
    
    % --- Algorithm 3: LMS ---
    y_lms = w_lms' * by;             % Filter Output
    e_lms = desired_sym - y_lms;     % Error
    w_lms = w_lms + mu_lms * conj(e_lms) * by; % Update
    
    xest_lms(n) = y_lms;
    err_lms_curve(n) = abs(e_lms)^2;
    
    % --- Algorithm 4: NLMS ---
    y_nlms = w_nlms' * by;           % Filter Output
    e_nlms = desired_sym - y_nlms;   % Error
    pwr = (by' * by) + 1e-5;         % Signal Power + epsilon
    w_nlms = w_nlms + (mu_nlms / pwr) * conj(e_nlms) * by; % Update
    
    xest_nlms(n) = y_nlms;
    err_nlms_curve(n) = abs(e_nlms)^2;
    
    % --- Algorithm 5: RLS (The New Addition) ---
    % A. Calculate Output
    y_rls = w_rls' * by;
    e_rls = desired_sym - y_rls;
    
    % B. Calculate Gain Vector (k)
    % Pi = P * u(n)
    Pi = P * by;
    k = Pi / (lambda_rls + by' * Pi);
    
    % C. Update Weights
    w_rls = w_rls + k * e_rls;
    
    % D. Update Inverse Correlation Matrix P
    P = (P - k * by' * P) / lambda_rls;
    
    xest_rls(n) = y_rls;
    err_rls_curve(n) = abs(e_rls)^2;
end

% =================================================================
% PART C: Results & Visualization
% =================================================================

% Plot Constellations (Converged - last 1000 symbols)
figure(2);
range = N-1000:N;

subplot(2,3,1); plot(xest_zf(range), '.'); title('Block: LS (ZF)'); axis([-2 2 -2 2]); axis square; grid on;
subplot(2,3,2); plot(xest_reg(range),'.'); title('Block: Robust LS'); axis([-2 2 -2 2]); axis square; grid on;
subplot(2,3,3); plot(y(range), 'r.');      title('Unequalized'); axis([-2 2 -2 2]); axis square; grid on;

subplot(2,3,4); plot(xest_lms(range), '.'); title('Adaptive: LMS'); axis([-2 2 -2 2]); axis square; grid on;
subplot(2,3,5); plot(xest_nlms(range),'.'); title('Adaptive: NLMS'); axis([-2 2 -2 2]); axis square; grid on;
subplot(2,3,6); plot(xest_rls(range), '.'); title('Adaptive: RLS'); axis([-2 2 -2 2]); axis square; grid on;

% Plot Convergence Speed (Learning Curve)
% We smooth the error curve to make it readable
window_smooth = 50;
figure(3);
semilogy(filter(ones(1,window_smooth)/window_smooth, 1, err_lms_curve), 'r', 'LineWidth', 1); hold on;
semilogy(filter(ones(1,window_smooth)/window_smooth, 1, err_nlms_curve), 'b', 'LineWidth', 1);
semilogy(filter(ones(1,window_smooth)/window_smooth, 1, err_rls_curve), 'k', 'LineWidth', 2);
legend('LMS', 'NLMS', 'RLS (Fastest)');
grid on;
xlabel('Symbol Index'); ylabel('Mean Squared Error (dB)');
title('Convergence Speed Comparison');
axis([0 1000 10^-3 10^1]);

% Final Console Output
fprintf('------------------------------------------\n');
fprintf('Method comparison complete.\n');
fprintf('Observe Figure 2 for constellation quality.\n');
fprintf('Observe Figure 3 for convergence speed.\n');
fprintf('Notice how RLS (Black line) drops to low error much faster than LMS.\n');