% === Updated RadioFM_stereo.m ===
clear all; close all; clc;

%% User parameters
centerFreq = 96.0e6; % RMF FM
basebandFs = 200e3;
plutoGain  = 70;
frameLen   = 4096*2;
audioFs    = 48000;
doStereo   = true;           % allow stereo decode (per-frame detection will enable/disable)
saveAudio  = true;           
audioFilename = 'pluto_fm_output.wav';
saveIQ     = false;          
iqFilename = 'pluto_fm_iq.dat';
deemphasisTau = 50e-6;       % Europe 50e-6

%% Create Pluto receiver System object
rx = comm.SDRRxPluto( ...
    'CenterFrequency', centerFreq, ...
    'BasebandSampleRate', basebandFs, ...
    'OutputDataType','double', ...
    'RadioID','usb:0', ...
    'GainSource','Manual', ...
    'Gain', plutoGain, ...
    'SamplesPerFrame', frameLen);

deviceWriter = audioDeviceWriter('SampleRate', audioFs);

%% Design filters
% Mono lowpass (0..15 kHz)
monoLP = designfilt('lowpassiir', 'PassbandFrequency', 15e3, ...
    'StopbandFrequency', 18e3, 'SampleRate', basebandFs, ...
    'PassbandRipple', 1, 'StopbandAttenuation', 60);

% Stereo bandpass to isolate L-R DSB around 38 kHz (keep 23..53 kHz)
stereoBP = designfilt('bandpassiir', 'FilterOrder', 8, ...
    'HalfPowerFrequency1', 23e3, 'HalfPowerFrequency2', 53e3, ...
    'SampleRate', basebandFs);

% Narrowband filter to extract 19 kHz pilot (±200 Hz or so)
pilotBP = designfilt('bandpassiir', 'FilterOrder', 6, ...
    'HalfPowerFrequency1', 18.5e3, 'HalfPowerFrequency2', 19.5e3, ...
    'SampleRate', basebandFs);

% De-emphasis coefficients (applied at audioFs)
alpha = (1/audioFs) / (deemphasisTau + 1/audioFs); % as in your original derivation
b_deem = alpha;
a_deem = [1 (alpha-1)];   % transfer function for y[n] = alpha*x[n] + (1-alpha)*y[n-1]

%% Misc storage
audioBuf = [];
if saveIQ
    fidIQ = fopen(iqFilename,'wb');
end

disp('Starting FM reception. Press Ctrl-C to stop.');

% runtime continuity variables
phaseOffset38 = 0;    % if using an explicit free-running carrier (as fallback)
prevSample = [];      % for conj multiply continuity if frames might be separated

% For pilot lock detection
pilotThreshold = 0.02; % empirical; might need tuning depending on signal strength

% fixed recording length (keeps your original behavior)
durationSec = 45;
tStart = tic;

try
    while toc(tStart) < durationSec
        rxIQ = rx();
        if isempty(rxIQ)
            pause(0.01);
            continue;
        end

        if saveIQ
            fwrite(fidIQ, [real(rxIQ) imag(rxIQ)] , 'double');
        end

        % ---------------------
        % FM demod (robust): angle of conj product
        % ---------------------
        % if we have previous sample from last frame, prepend for continuity
        if ~isempty(prevSample)
            rx_cont = [prevSample; rxIQ];
            d = rx_cont(2:end) .* conj(rx_cont(1:end-1));
            prevSample = rxIQ(end); % update for next frame continuity
            % drop first demod sample (corresponds to prevSample->first sample)
            d = d(2:end);
        else
            d = rxIQ(2:end) .* conj(rxIQ(1:end-1));
            prevSample = rxIQ(end);
        end
        dAng = angle(d);                 % radians per sample
        % convert to normalized audio baseband. scale by basebandFs/(2*pi*maxDeviation)
        maxDeviation = 75e3;             % typical WFM deviation
        instFreqHz = (basebandFs/(2*pi)) * dAng;  % Hz (instantaneous deviation)
        audioBase = instFreqHz / maxDeviation;    % normalized approximately within ±1

        % Ensure audioBase length matches frameLen (we lost one sample in diff)
        % pad with the first value to restore frame length
        audioBase = [audioBase(1); audioBase];

        % ---------------------
        % Mono (L+R)
        % ---------------------
        mono = filtfilt(monoLP, audioBase);  % zero-phase mono

        % ---------------------
        % Pilot detection (19 kHz)
        % ---------------------
        pilot = filtfilt(pilotBP, audioBase); % extract pilot tone
        % get analytic signal to measure phase and amplitude
        pilotAnalytic = hilbert(pilot);
        pilotAmp = abs(pilotAnalytic);
        pilotPhase = angle(pilotAnalytic);   % per-sample instantaneous phase ~= 2π*19k*t + phi

        % measure pilot strength (RMS or mean envelope)
        pilotStrength = mean(pilotAmp);

        if doStereo && (pilotStrength > pilotThreshold)
            % ---- stereo path ----
            % isolate stereo subband (23..53k)
            stereoSub = filtfilt(stereoBP, audioBase);

            % Use pilot phase to create coherent 38 kHz carrier:
            % pilotPhase ~ 2π*19k*t + φ  -> doubling gives ~ 2π*38k*t + 2φ
            % create carrier = 2*cos(2*pilotPhase). Doing this per-sample keeps phase lock.
            % Note: pilotPhase may have noise; smoothing or PLL can improve results.
            carrier38 = 2 .* cos(2 .* pilotPhase);   % length = frameLen

            % multiply (demodulate DSB-SC L-R to baseband)
            demodLR = stereoSub .* carrier38;

            % lowpass to extract L-R (audio band)
            LR = filtfilt(monoLP, demodLR);

            % left and right
            left  = 0.5 * (mono + LR);
            right = 0.5 * (mono - LR);

            % resample each channel to audioFs
            leftRes  = resample(left,  audioFs, basebandFs);
            rightRes = resample(right, audioFs, basebandFs);
            audioFrame = [leftRes, rightRes];

        else
            % ---- mono fallback ----
            monoRes = resample(mono, audioFs, basebandFs);
            audioFrame = [monoRes, monoRes]; % duplicate for stereo playback
        end

        % De-emphasis (IIR) applied across channels
        audioFrame(:,1) = filter(b_deem, a_deem, audioFrame(:,1));
        audioFrame(:,2) = filter(b_deem, a_deem, audioFrame(:,2));

        % Play
        deviceWriter(audioFrame);

        % Append
        if saveAudio
            audioBuf = [audioBuf; audioFrame]; %#ok<AGROW>
        end
    end

catch ME
    release(rx);
    release(deviceWriter);
    if saveIQ
        fclose(fidIQ);
    end
    warning('Reception stopped: %s', ME.message);
end

if saveAudio && ~isempty(audioBuf)
    audioBuf = audioBuf / max(abs(audioBuf(:)));
    audiowrite(audioFilename, audioBuf, audioFs);
    fprintf('Audio saved to %s\n', audioFilename);
end

disp('Done.');
