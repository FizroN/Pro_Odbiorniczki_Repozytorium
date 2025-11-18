%% Exercise 20.7 – Using pulse shaping filters (RC & SRRC)
clear; close all; clc;

%% ==================== PARAMETERS ===========================
K = 6;          % samples per symbol (interpolation factor)
Nsym = 300;     % QAM symbol count
r = 0.35;       % rolloff
Ks = 6;         % number of symbols in pulse shaping filter
Npsf = K * Ks;  % filter length
if mod(Npsf,2) ~= 0
    Npsf = Npsf + 1;
end
n = (-Npsf/2 : Npsf/2);     % symmetric index

%% ==================== GENERATE SYMBOL STREAM =================
Ik = sign(randn(1,Nsym));
Qk = sign(randn(1,Nsym));
IQk = Ik + 1j*Qk;

%% ==================== RC FILTER (analytical) =================
num = sin(pi*n/K) ./ (pi*n/K);
num(isnan(num)) = 1;   % sinc(0)=1

den = 1 - (2*r*n/K).^2;
rc = num .* cos(pi*r*n/K) ./ den;
% singularities
idx = find(abs(den) < 1e-12);
for i=idx
    rc(i) = r/2 * sin(pi/(2*r));
end
hRC = rc / sum(rc);     % normalize

%% ==================== SRRC FILTER ============================
% standard analytical formula
x = n/K;
hSRRC = (4*r/pi) * (cos((1+r)*pi*x) + (sin((1-r)*pi*x) ./ (4*r*x))) ./ (1 - (4*r*x).^2);
hSRRC(n==0) = (1 + r*(4/pi - 1));        % center sample
hSRRC = hSRRC / sum(hSRRC);              % normalize

%% ==================== INTERPOLATION ==========================
% Zero-stuffing
IQ0 = zeros(1, length(IQk)*K);
IQ0(1:K:end) = IQk;

% Interpolation by convolution
IQ_RC   = conv(IQ0, hRC,   'same');
IQ_SRRC = conv(IQ0, hSRRC, 'same');

%% ==================== TIME DOMAIN PLOTS ======================
figure;
subplot(2,1,1);
stem(real(IQ_RC(1:250)), 'filled');
title('I(n) interpolated with RC'); grid on;

subplot(2,1,2);
stem(imag(IQ_RC(1:250)), 'filled');
title('Q(n) interpolated with RC'); grid on;

figure;
subplot(2,1,1);
stem(real(IQ_SRRC(1:250)), 'filled');
title('I(n) interpolated with SRRC'); grid on;

subplot(2,1,2);
stem(imag(IQ_SRRC(1:250)), 'filled');
title('Q(n) interpolated with SRRC'); grid on;

%% ==================== PHASOR DIAGRAMS ========================
figure;
plot(real(IQ_RC), imag(IQ_RC), '.');
title('Phasor diagram – RC'); axis equal; grid on;

figure;
plot(real(IQ_SRRC), imag(IQ_SRRC), '.');
title('Phasor diagram – SRRC'); axis equal; grid on;

%% ==================== EYE DIAGRAM ============================
L = K;                  % symbol duration in samples
offset = 3*K;           % skip filter transient
sig = IQ_SRRC(offset:end);

% make number of samples divisible by L
M = floor(length(sig)/L)*L;
sig = sig(1:M);

E = reshape(sig, L, []);

figure;
plot(real(E), 'LineWidth', 1);
title('Eye diagram – SRRC (I component)');
xlabel('Sample'); grid on;

figure;
plot(imag(E), 'LineWidth', 1);
title('Eye diagram – SRRC (Q component)');
xlabel('Sample'); grid on;
