% pluto_tx.m
% Transmits a signal using the ADALM-PLUTO SDR.
% Based on lab21_ex_receiver.m

clear all; close all;

%% 1. Parameters (MUST MATCH THE RECEIVER)
% --- SDR Parameters ---
f_carrier = 915e6;         % 915 MHz (Check local regulations!)
pluto_ip = 'ip:192.168.2.1'; % !!! CHANGE THIS to your Pluto's IP

% --- Waveform Parameters ---
modtype = 'QPSK';      % QPSK is robust. 4QAM, 16QAM also work.
Rs = 57e3;             % Symbol Rate (57 kS/s)
K = 8;                 % Samples per Symbol
Ns = 8;
fs = Rs * K;           % Baseband Sample Rate (456 kS/s)
r = 0.35;              % Roll-off factor for pulse shaping
rctype = 'sqrt';       % Root-Raised Cosine
Ndata = 1000;          % Number of data symbols

%% 2. Generate the "Frame" to Transmit (from lab21)
% This section creates the data payload and the sync header.
Npsf = K*Ns+1; Mpsf = (Npsf-1)/2;
[IQcodes, Nstates, Nbits, R ] = IQdef( modtype );

% IQk of Header (The known sequence the receiver will look for)
[numHead, Nhead ] = modtype2header( modtype );
IQkHead = numbers2IQ( numHead, modtype, IQcodes );

% IQk of Data (The random data we want to send)
numData = floor( Nstates*(rand(Ndata,1)-10*eps) );
IQkData = numbers2IQ( numData, modtype, IQcodes );

% --- Build the full frame ---
% We'll send: [Header | Data]
% (We can add a "tail" of zeros for the matched filter to flush)
IQk = [ IQkHead  IQkData ];
num = [ numHead' numData' ];

%% 3. Pulse Shaping (from lab21)
% Turn the digital symbols (IQk) into a smooth, band-limited
% complex baseband signal (IQn).
[IQn, hpsf ] = IQ2psf( IQk, K, Ns, r, rctype );

% --- Prepare signal for SDR ---
% 1. Scale to [-1, 1] to prevent clipping
max_amp = max(abs(IQn));
tx_signal = IQn / max_amp;

% 2. Must be a column vector for the SDR object
tx_signal = tx_signal.'; 

% Optional: Plot the constellation of the signal we're about to send
N = length( IQn ); n = Npsf : N-Npsf+1; ns = Npsf : K : N-Npsf+1;
figure; 
plot( real(IQn(n)), imag(IQn(n)), real(IQn(ns)), imag(IQn(ns)),'ro','MarkerFaceColor','red'); 
grid; title('TX Signal Constellation (Baseband)');
xlabel('I'); ylabel('Q');

%% 4. Set up ADALM-PLUTO Transmitter
disp('Setting up Pluto TX...');
try
    tx = sdrtx('Pluto', ...
        'RadioID', pluto_ip, ...
        'CenterFrequency', f_carrier, ...
        'BasebandSampleRate', fs, ...
        'Gain', -10); % Start with low gain (-10 dB) to avoid issues

    % --- Transmit! ---
    % We use transmitRepeat to send the frame over and over.
    % This makes it much easier for the receiver to lock on.
    disp('Transmitting... Press Ctrl+C to stop.');
    transmitRepeat(tx, tx_signal);

catch ME
    disp('Error setting up or running Pluto TX:');
    disp(ME.message);
    if (exist('tx', 'var'))
        release(tx);
    end
end

% Note: The script will loop on transmitRepeat() until you
% press Ctrl+C. You must stop this script before running the receiver
% if you only have one Pluto, or run them on two separate computers.
% If you have two Plutos, you can run them at the same time.
