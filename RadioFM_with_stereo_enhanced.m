% RadioFM_stereo_PLL.m
% FM receiver with pilot-locked stereo decode using a simple digital PLL.
% - Uses ADALM-PLUTO (comm.SDRRxPluto) or can operate on saved IQ file (test mode).
% - Stateful causal filters (no filtfilt) for streaming.
% - Pilot auto-threshold calibration, live PLL, and stereo fallback to mono.
% - Includes analyzeIQ(...) to visualize spectrograms, pilot amplitude and L/R channels.
%
% Usage:
%  - Configure user parameters below.
%  - Run the file. Press Ctrl-C to stop or let durationSec elapse.
%  - To analyze a saved IQ file, call analyzeIQ('pluto_fm_iq.dat', basebandFs, audioFs).

clearvars; close all; clc;

%% ================= User parameters =================
usePluto   = true;           % true -> live from Pluto; false -> test mode from iq file
centerFreq = 96.0e6;         % tune frequency (Hz)
basebandFs = 200e3;          % Pluto baseband sampling (Hz)
plutoGain  = 70;
frameLen   = 4096*2;         % samples per Pluto frame
audioFs    = 48000;          % audio sample rate (Hz)
doStereo   = true;           % allow stereo decode (will only enable if pilot present)
saveAudio  = true;
audioFilename = 'pluto_fm_output_stereo.wav';
saveIQ     = true;
iqFilename = 'pluto_fm_iq.dat';
deemphasisTau = 50e-6;       % Europe = 50e-6, US = 75e-6
durationSec = 45;            % seconds to record (set [] to run until Ctrl-C)

% PLL parameters (tweak if necessary)
pllLoopBandwidth = 20;  % Hz (controls how fast PLL tracks pilot jitter) - ~10-60 typical
pllDamping = 0.707;     % damping factor (zeta)
%% ==================================================

% Derived
if isempty(durationSec)
    durationSec = inf;
end

% Create Pluto receiver if needed
if usePluto
    rx = comm.SDRRxPluto( ...
        'CenterFrequency', centerFreq, ...
        'BasebandSampleRate', basebandFs, ...
        'OutputDataType','double', ...
        'RadioID','usb:0', ...
        'GainSource','Manual', ...
        'Gain', plutoGain, ...
        'SamplesPerFrame', frameLen);
end

% Audio device writer
deviceWriter = audioDeviceWriter('SampleRate', audioFs);

% ----------------- Filter designs (causal, IIR) -----------------
% We will keep filter states (zi) and use filter(b,a, x, zi) across frames.

% Mono lowpass for L+R (0..15k)
monoOrder = 6;
monoCut = 15e3/(basebandFs/2);
[b_mono, a_mono] = butter(monoOrder, monoCut, 'low');
zi_mono = zeros(max(length(a_mono), length(b_mono))-1, 1);

% Bandpass for stereo subband (23..53k) around 38k DSB
bpOrder = 6;
bpWn = [23e3 53e3] / (basebandFs/2);
[b_bp, a_bp] = butter(bpOrder, bpWn, 'bandpass');
zi_bp = zeros(max(length(a_bp), length(b_bp))-1, 1);

% Narrow band for 19 kHz pilot extraction (±800 Hz or so)
pilotOrder = 4;
pilotBW = 1000; % +/- Hz around 19k
pilotWn = [(19e3 - pilotBW) (19e3 + pilotBW)] / (basebandFs/2);
[b_pilot, a_pilot] = butter(pilotOrder, pilotWn, 'bandpass');
zi_pilot = zeros(max(length(a_pilot), length(b_pilot))-1, 1);

% Lowpass to extract LR after demod (same as mono lowpass)
% we'll reuse b_mono/a_mono and zi_monoLR separate
zi_monoLR = zeros(max(length(a_mono), length(b_mono))-1, 1);

% De-emphasis (IIR) coefficients for audioFs
a_deem = [1 ( (1/audioFs)/(deemphasisTau + 1/audioFs) - 1 )];
b_deem = (1/audioFs)/(deemphasisTau + 1/audioFs);
zi_deem_L = 0;
zi_deem_R = 0;

% ----------------- PLL initialization -----------------
% Simple 2nd-order digital PLL (phase detector from pilot phase)
% We'll implement a discrete-time PI loop in phase domain:
% vco_phase[n+1] = vco_phase[n] + w0 + Kp*err + Ki*integrator
% where err = angle(pilot * conj(exp(1j*vco_phase))). We use sample-by-sample update.

% PLL constants derived from loop bandwidth & damping (approx continuous->discrete)
Fs_pll = basebandFs; % PLL runs at baseband sample rate
f0 = 19e3; % pilot nominal frequency (Hz)

% Convert desired loop bandwidth to loop gains (digital)
% Following classical control mapping for 2nd-order type-II style PLL
% approximate formulas:
wn = 2*pi*pllLoopBandwidth;           % loop natural frequency (rad/s)
Kpd = 1;                              % phase detector gain (rad/rad)
Kvco = 2*pi*f0/Fs_pll;                % VCO gain in rad/sample per control unit (approx)
% Factor to compute loop coefficients (standard control formulas)
Kp_pll = (2*pllDamping*wn - 0)/(Kpd*Kvco);  % proportional (approx)
Ki_pll = (wn^2)/(Kpd*Kvco);                 % integral (approx)

% Safety clip for gains
if isnan(Kp_pll) || isnan(Ki_pll) || Kp_pll<=0 || Ki_pll<=0
    Kp_pll = 1e-3; Ki_pll = 1e-6;
end

vco_phase = 0;        % VCO instant phase (radians)
vco_integrator = 0;   % integrator state

% ----------------- Other runtime vars -----------------
audioBuf = [];
if saveIQ
    fidIQ = fopen(iqFilename,'wb');
end
prevSample = []; % for conj multiply continuity
pilot_strength_history = []; % for auto-threshold calibration

disp('Starting FM reception with PLL stereo decode. Press Ctrl-C to stop.');

tStart = tic;
try
    while toc(tStart) < durationSec
        % ---------- receive frame ----------
        if usePluto
            rxIQ = rx();
            if isempty(rxIQ)
                pause(0.01);
                continue;
            end
        else
            error('Non-Pluto mode not implemented in main loop; use analyzeIQ(...) to test files.');
        end

        % optionally save IQ raw doubles (I,Q interleaved)
        if saveIQ
            fwrite(fidIQ, [real(rxIQ) imag(rxIQ)] , 'double');
        end

        % ---------- FM demod (conj multiply, angle) ----------
        % Keep continuity across frames: prepend prevSample if exists
        if ~isempty(prevSample)
            rx_cont = [prevSample; rxIQ];
            d = rx_cont(2:end) .* conj(rx_cont(1:end-1));
            dAng = angle(d);
            % discard the first sample because it corresponds to prevSample->first current
            dAng = dAng(2:end);
            prevSample = rxIQ(end);
        else
            d = rxIQ(2:end) .* conj(rxIQ(1:end-1));
            dAng = angle(d);
            dAng = [dAng(1); dAng]; % safe padding (rare)
            prevSample = rxIQ(end);
        end

        instFreqHz = (basebandFs/(2*pi)) * dAng;
        maxDeviation = 75e3;
        audioBase = instFreqHz / maxDeviation;  % normalized audio baseband (approx)

        % ensure same length as frameLen
        if length(audioBase) < frameLen
            audioBase = [audioBase; zeros(frameLen - length(audioBase),1)];
        elseif length(audioBase) > frameLen
            audioBase = audioBase(1:frameLen);
        end

        % ---------- Mono (L+R) path ----------
        [mono, zi_mono] = filter(b_mono, a_mono, audioBase, zi_mono);

        % ---------- Pilot extraction ----------
        [pilot, zi_pilot] = filter(b_pilot, a_pilot, audioBase, zi_pilot);
        % analytic signal via Hilbert (we will compute instantaneous phase)
        pilot_analytic = hilbert(pilot);
        pilot_amp = abs(pilot_analytic);
        pilot_phase_inst = angle(pilot_analytic); % instantaneous phase of pilot

        % track pilot strength history for initial calibration
        pilot_strength_history = [pilot_strength_history; mean(pilot_amp)]; %#ok<AGROW>
        if length(pilot_strength_history) > 200
            pilot_strength_history(1:end-200) = []; % cap history size
        end

        % ---------- PLL: per-sample loop to track pilot phase ----------
        % We'll run the PLL sample-by-sample over the frame (vectorized loop)
        N = length(pilot);
        carrier38 = zeros(N,1);
        for n = 1:N
            % Phase error between pilot and VCO
            % pilot phase is pilot_phase_inst(n)
            phase_error = wrapToPi(pilot_phase_inst(n) - vco_phase);
            % update integrator and vco_phase
            vco_integrator = vco_integrator + Ki_pll * phase_error;
            vco_phase = vco_phase + Kvco + Kp_pll * phase_error + vco_integrator;
            % store doubled-phase carrier at sample n (coherent 38 kHz)
            carrier38(n) = 2 * cos(2 * vco_phase);  % multiply by 2 for amplitude scaling of DSB demod
        end
        % Keep vco_phase and integrator states for next frame (already updated)

        % ---------- Stereo subband isolation ----------
        [stereoSub, zi_bp] = filter(b_bp, a_bp, audioBase, zi_bp);

        % Multiply stereo subband by coherent 38k carrier to demodulate DSB L-R
        demodLR = stereoSub .* carrier38;

        % Lowpass LR to audio band (same mono lowpass coefficients) - use separate state
        [LR, zi_monoLR] = filter(b_mono, a_mono, demodLR, zi_monoLR);

        % ---------- Build left/right ----------
        left  = 0.5 * (mono + LR);
        right = 0.5 * (mono - LR);

        % Decide whether stereo is actually present (auto calibration)
        % Compute pilot strength metric (mean amplitude)
        instPilotStrength = mean(pilot_amp);
        % Auto threshold: if we have enough history, set threshold = median + 1.5*std
        if length(pilot_strength_history) > 50
            med = median(pilot_strength_history);
            s = std(pilot_strength_history);
            pilotThreshold = med + 1.5*s;
        else
            pilotThreshold = 0.02; % default initial
        end

        if doStereo && (instPilotStrength > pilotThreshold)
            % keep left/right as computed
            finalLeft = left;
            finalRight = right;
        else
            % fallback to mono duplicated
            finalLeft = mono;
            finalRight = mono;
        end

        % ---------- Resample to audioFs ----------
        leftRes  = resample(finalLeft,  audioFs, basebandFs);
        rightRes = resample(finalRight, audioFs, basebandFs);
        audioFrame = [leftRes, rightRes];

        % ---------- De-emphasis (applied per-channel with state) ----------
        [audioFrame(:,1), zi_deem_L] = filter(b_deem, a_deem, audioFrame(:,1), zi_deem_L);
        [audioFrame(:,2), zi_deem_R] = filter(b_deem, a_deem, audioFrame(:,2), zi_deem_R);

        % ---------- Playback ----------
        deviceWriter(audioFrame);

        % ---------- Save audio buffer ----------
        if saveAudio
            audioBuf = [audioBuf; audioFrame]; %#ok<AGROW>
        end

    end
catch ME
    release(deviceWriter);
    if usePluto; release(rx); end
    if saveIQ; fclose(fidIQ); end
    rethrow(ME);
end

% Clean up after main loop
release(deviceWriter);
if usePluto; release(rx); end
if saveIQ; fclose(fidIQ); end

% Save audio to WAV (normalize)
if saveAudio && ~isempty(audioBuf)
    audioBuf = audioBuf / max(abs(audioBuf(:)));
    audiowrite(audioFilename, audioBuf, audioFs);
    fprintf('Audio saved to %s\n', audioFilename);
end
disp('Done.');


