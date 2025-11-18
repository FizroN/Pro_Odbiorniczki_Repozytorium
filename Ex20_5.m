%% Exercise 20.5 — RC filter as convolution of two SRRC filters (firrcos version)
clear; close all; clc;

r      = 0.35;      % rolloff
Nsymb  = 24;        % span in symbols
Ksymb  = 8;         % samples per symbol
Fs     = 1;         % sample time = 1

%% ----------------------------------------------------
%         RC and SRRC design using firrcos()
% ----------------------------------------------------
try
    % RC — raised cosine
    hRC = firrcos(Nsymb*Ksymb, ...          % filter order
                  1/(2*Nsymb), ...          % cutoff freq (symbol rate/2)
                  r, ...                    % rolloff
                  Fs, ...
                  'rolloff', 'normal');     % standard RC

    % SRRC — square-root raised cosine
    hSRRC = firrcos(Nsymb*Ksymb, ...
                    1/(2*Nsymb), ...
                    r, ...
                    Fs, ...
                    'type', 'sqrt');

catch
    warning('firrcos() not found — using rcosdesign() instead');
    % Equivalent design using Communications Toolbox
    hRC   = rcosdesign(r, Nsymb, Ksymb, 'normal');
    hSRRC = rcosdesign(r, Nsymb, Ksymb, 'sqrt');
end

%% Normalize both to unit energy (good practice)
hRC   = hRC(:)' / sum(hRC);
hSRRC = hSRRC(:)' / sum(hSRRC);

%% ----------------------------------------------------
%               Convolution: SRRC * SRRC
% ----------------------------------------------------
hConv = conv(hSRRC, hSRRC);

% Lengths
Lr   = length(hRC);
Lc   = length(hConv);

% Center indices
center_RC   = ceil(Lr/2);
center_Conv = ceil(Lc/2);

% Trim convolution to match RC length (center aligned)
start_idx   = center_Conv - (center_RC-1);
end_idx     = start_idx + Lr - 1;
hConv_trim  = hConv(start_idx:end_idx);

% Normalize energies to match
hConv_trim = hConv_trim * (sum(hRC)/sum(hConv_trim));

%% ----------------------------------------------------
%                  Time-domain comparison
% ----------------------------------------------------
figure;
subplot(2,1,1);
plot(hRC,'b','LineWidth',1.2); hold on;
plot(hConv_trim,'r--','LineWidth',1.2);
title('Time Domain: RC vs SRRC ⊗ SRRC');
legend('h_{RC}','conv(h_{SRRC},h_{SRRC})');
grid on;

%% ----------------------------------------------------
%                Frequency-domain comparison
% ----------------------------------------------------
Nfft = 16384;
f = linspace(-0.5,0.5,Nfft);

Hrc = fftshift(fft(hRC, Nfft));
Hc  = fftshift(fft(hConv_trim, Nfft));

subplot(2,1,2);
plot(f,20*log10(abs(Hrc)+eps),'b','LineWidth',1.2); hold on;
plot(f,20*log10(abs(Hc)+eps),'r--','LineWidth',1.2);
title('Frequency Domain: |H_{RC}| vs |H_{SRRC} * H_{SRRC}|');
xlabel('Normalized frequency'); ylabel('Magnitude (dB)');
grid on;
axis([-0.5 0.5 -80 5]);
legend('|H_{RC}|','|H_{SRRC⊗SRRC}|');

%% ----------------------------------------------------
%                     Numerical Error
% ----------------------------------------------------
err = max(abs(hRC - hConv_trim));
disp(['Max abs error between RC and SRRC-convolution = ' num2str(err)]);

