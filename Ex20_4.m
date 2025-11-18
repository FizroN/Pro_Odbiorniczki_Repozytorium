%% Exercise 20.4 — Testing RC/SRRC filter length
% r = 0.35; Nsymb = 24;
% Ksymb list = 2, 4, 6, 8, 10, 12
% One figure for RC, one for SRRC

clear; close all; clc;

% Given parameters
r = 0.35;
Nsymb = 24;
K_list = [2 4 6 8 10 12];
colors = lines(length(K_list));

%% ============================================================
%                 FIGURE 1 — RC FILTERS
% ============================================================

figure('Name','RC Filters for different Ksymb','NumberTitle','off');

% Impulse responses
subplot(2,1,1); hold on; grid on;
title('Impulse Responses — Raised Cosine');
xlabel('n - N_c'); ylabel('h_{RC}[n]');

% Magnitude responses
subplot(2,1,2); hold on; grid on;
title('Magnitude Responses |H_{RC}(f)| (dB)');
xlabel('Normalized Frequency'); ylabel('Magnitude (dB)');
axis([-0.5 0.5 -80 5]);

for k = 1:length(K_list)
    
    Ksymb = K_list(k);
    Npsf = Nsymb * Ksymb;
    n = 0:Npsf;
    Nc = Npsf/2;

    % ----------- RC pulse definition -----------
    h_rc = sin(pi*(n-Nc)/Nsymb) ./ (pi*(n-Nc)/Nsymb) .* ...
           cos(pi*r*(n-Nc)/Nsymb) ./ (1 - (2*r*(n-Nc)/Nsymb).^2);
    h_rc(Nc+1) = 1;

    n0 = find(abs(1-(2*r*(n-Nc)/Nsymb).^2) < sqrt(eps));
    h_rc(n0) = r/2 * sin(pi/(2*r));

    h_rc = h_rc / Nsymb;

    % ----------- Impulse plot -----------
    subplot(2,1,1);
    plot(n-Nc, h_rc, 'Color', colors(k,:), 'LineWidth', 1.3);

    % ----------- Magnitude response -----------
    subplot(2,1,2);
    Nfft = 32768;
    H = fftshift(fft(h_rc, Nfft));
    f = linspace(-0.5, 0.5, Nfft);
    plot(f, 20*log10(abs(H)+eps), 'Color', colors(k,:), 'LineWidth', 1.3);
end

subplot(2,1,1); legend("K=2","K=4","K=6","K=8","K=10","K=12");
subplot(2,1,2); legend("K=2","K=4","K=6","K=8","K=10","K=12");


%% ============================================================
%                 FIGURE 2 — SRRC FILTERS
% ============================================================

figure('Name','SRRC Filters for different Ksymb','NumberTitle','off');

% Impulse responses
subplot(2,1,1); hold on; grid on;
title('Impulse Responses — Square Root Raised Cosine');
xlabel('n - N_c'); ylabel('h_{SRRC}[n]');

% Magnitude responses
subplot(2,1,2); hold on; grid on;
title('Magnitude Responses |H_{SRRC}(f)| (dB)');
xlabel('Normalized Frequency'); ylabel('Magnitude (dB)');
axis([-0.5 0.5 -80 5]);

for k = 1:length(K_list)
    
    Ksymb = K_list(k);
    Npsf = Nsymb * Ksymb;
    n = 0:Npsf;
    Nc = Npsf/2;

    % ------------- SRRC pulse definition -------------
    n0 = find(abs(1-(4*r*(n-Nc)/Nsymb).^2) < 5*eps);
    h_srrc = 4*r/pi * ( ...
        cos(pi*(1+r)*(n-Nc)/Nsymb) + ...
        (Nsymb./(4*r*(n-Nc))) .* sin(pi*(1-r)*(n-Nc)/Nsymb) ...
        ) ./ (1 - (4*r*(n-Nc)/Nsymb).^2);

    h_srrc(Nc+1) = (1 + r*(4/pi - 1));
    h_srrc = h_srrc / Nsymb;

    % ------------- Impulse plot -------------
    subplot(2,1,1);
    plot(n-Nc, h_srrc, 'Color', colors(k,:), 'LineWidth', 1.3);

    % ------------- Magnitude plot -------------
    subplot(2,1,2);
    Nfft = 32768;
    H = fftshift(fft(h_srrc, Nfft));
    f   = linspace(-0.5, 0.5, Nfft);
    plot(f, 20*log10(abs(H)+eps), 'Color', colors(k,:), 'LineWidth', 1.3);
end

subplot(2,1,1); legend("K=2","K=4","K=6","K=8","K=10","K=12");
subplot(2,1,2); legend("K=2","K=4","K=6","K=8","K=10","K=12");
