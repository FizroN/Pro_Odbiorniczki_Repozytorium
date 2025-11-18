%% Exercise 20.7 - Using pulse shaping filters (RC & SRRC)
clear; close all; clc;

%% ------------------- PARAMETERS ------------------------
K = 6;           % interpolation factor (samples per symbol)
Nsym = 300;      % number of QAM symbols

%% ------------------- GENERATE SYMBOL STREAM (program 20.2 style) ----
% 4-QAM (QPSK)
Ik = sign(randn(1, Nsym));
Qk = sign(randn(1, Nsym));
IQk = Ik + 1j*Qk;

%% ------------------- LOAD PREVIOUS FILTERS --------------------------
load hRC.mat;      % <- save from Exercise 20.6
load hSRRC.mat;

%% ------------------- INTERPOLATION METHOD 1 (upfirdown) -------------
IQn_RC   = upfirdown(IQk,   hRC,   K, 1);
IQn_SRRC = upfirdown(IQk, hSRRC, K, 1);

%% ------------------- INTERPOLATION METHOD 2 (manual) ----------------
% Zero insertion
IQ0 = zeros(1, Nsym*K);
IQ0(1:K:end) = IQk;

% Filter
IQn2_RC   = conv(IQ0,   hRC,   'same');
IQn2_SRRC = conv(IQ0, hSRRC, 'same');

%% ------------------- PLOTS: TIME DOMAIN ------------------------------
figure;
subplot(2,1,1);
stem(real(IQn_RC(1:200)), 'filled');
title('I(n) interpolated with RC'); grid on;

subplot(2,1,2);
stem(imag(IQn_RC(1:200)), 'filled');
title('Q(n) interpolated with RC'); grid on;

figure;
subplot(2,1,1);
stem(real(IQn_SRRC(1:200)), 'filled');
title('I(n) interpolated with SRRC'); grid on;

subplot(2,1,2);
stem(imag(IQn_SRRC(1:200)), 'filled');
title('Q(n) interpolated with SRRC'); grid on;

%% ------------------- PHASOR DIAGRAMS --------------------------------
figure;
plot(real(IQn_RC), imag(IQn_RC), '.');
title('Phasor Diagram (RC)'); xlabel('I'); ylabel('Q'); grid on; axis equal;

figure;
plot(real(IQn_SRRC), imag(IQn_SRRC), '.');
title('Phasor Diagram (SRRC)'); xlabel('I'); ylabel('Q'); grid on; axis equal;

%% ------------------- OPTIONAL: EYE DIAGRAM --------------------------
% Extract only samples aligned to the symbol period
L = K;   % one symbol = K samples

% Remove filter transient -- skip first 3L samples
offset = 3*L;
seq = IQn_RC(offset:end);

% Truncate to full eye diagram matrix
N = floor(length(seq)/L) * L;
seq = seq(1:N);

% Reshape: each row = one eye
eyeM = reshape(seq, L, []);

figure;
plot(real(eyeM));
title('Eye Diagram (RC)'); xlabel('Sample'); ylabel('I'); grid on;

figure;
plot(imag(eyeM));
title('Eye Diagram (RC)'); xlabel('Sample'); ylabel('Q'); grid on;
