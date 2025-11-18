%% lab20_ex_rcos_srrc_separate.m
% RC and SRRC for r = 0.1, 0.5, 0.9 — SEPARATE FIGURES

clear; close all; clc;

Nsymb = 10;      % samples per symbol
Ksymb = 4;       % symbols spanned by filter
Npsf  = Nsymb*Ksymb;

n  = 0:Npsf;
Nc = Npsf/2;

rolloff_list = [0.1 0.5 0.9];
colors = ['r','g','b'];

%% ==============================================================
%                     FIGURE 1 — RC ONLY
% ==============================================================

figure('Name','Raised Cosine (RC) — r = 0.1, 0.5, 0.9','NumberTitle','off');

% Impulse response subplot
subplot(2,1,1); hold on; grid on;
title('RC Impulse Responses');
xlabel('n - N_c'); ylabel('h_{RC}[n]');

% Magnitude response subplot
subplot(2,1,2); hold on; grid on;
title('RC Magnitude Responses |H(f)| (dB)');
xlabel('Normalized Frequency'); ylabel('Magnitude (dB)');
axis([-0.5 0.5 -80 5]);

for k = 1:length(rolloff_list)
    r = rolloff_list(k);

    % RC pulse
    h_rc = sin(pi*(n-Nc)/Nsymb) ./ (pi*(n-Nc)/Nsymb) .* ...
           cos(pi*r*(n-Nc)/Nsymb) ./ (1 - (2*r*(n-Nc)/Nsymb).^2);
    h_rc(Nc+1) = 1;
    n0 = find(abs(1-(2*r*(n-Nc)/Nsymb).^2) < sqrt(eps));
    h_rc(n0) = r/2 * sin(pi/(2*r));
    h_rc = h_rc / Nsymb;

    % Plot impulse response
    subplot(2,1,1);
    plot(n-Nc, h_rc, 'Color', colors(k), 'LineWidth', 1.3);

    % Plot magnitude
    subplot(2,1,2);
    Nfft = 16384;
    Hrc = fftshift(fft(h_rc, Nfft));
    f = linspace(-0.5, 0.5, Nfft);
    plot(f, 20*log10(abs(Hrc)+eps), 'Color', colors(k), 'LineWidth', 1.3);
end

subplot(2,1,1);
legend('r=0.1','r=0.5','r=0.9');

subplot(2,1,2);
legend('r=0.1','r=0.5','r=0.9');


%% ==============================================================
%                     FIGURE 2 — SRRC ONLY
% ==============================================================

figure('Name','Square-Root RC (SRRC) — r = 0.1, 0.5, 0.9','NumberTitle','off');

% Impulse response subplot
subplot(2,1,1); hold on; grid on;
title('SRRC Impulse Responses');
xlabel('n - N_c'); ylabel('h_{SRRC}[n]');

% Magnitude response subplot
subplot(2,1,2); hold on; grid on;
title('SRRC Magnitude Responses |H(f)| (dB)');
xlabel('Normalized Frequency'); ylabel('Magnitude (dB)');
axis([-0.5 0.5 -80 5]);

for k = 1:length(rolloff_list)
    r = rolloff_list(k);

    % SRRC pulse
    n0 = find(abs(1-(4*r*(n-Nc)/Nsymb).^2) < 5*eps);
    if ~isempty(n0)
        warning('SRRC singularity near pulse center.');
    end
    
    h_srrc = 4*r/pi * ( ...
        cos(pi*(1+r)*(n-Nc)/Nsymb) + ...
        (Nsymb./(4*r*(n-Nc))) .* sin(pi*(1-r)*(n-Nc)/Nsymb) ...
        ) ./ (1 - (4*r*(n-Nc)/Nsymb).^2);

    h_srrc(Nc+1) = (1 + r*(4/pi - 1));
    h_srrc = h_srrc / Nsymb;

    % Impulse response
    subplot(2,1,1);
    plot(n-Nc, h_srrc, 'Color', colors(k), 'LineWidth', 1.3);

    % Magnitude
    subplot(2,1,2);
    Nfft = 16384;
    Hsrrc = fftshift(fft(h_srrc, Nfft));
    f = linspace(-0.5, 0.5, Nfft);
    plot(f, 20*log10(abs(Hsrrc)+eps), 'Color', colors(k), 'LineWidth', 1.3);
end

subplot(2,1,1);
legend('r=0.1','r=0.5','r=0.9');

subplot(2,1,2);
legend('r=0.1','r=0.5','r=0.9');
