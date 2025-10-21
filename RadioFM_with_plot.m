% RadioFM.m
% FM receiver using ADALM-PLUTO and MATLAB
% - Mono by default, optional stereo decode (set doStereo = true)
% - Saves audio to WAV and optionally IQ to file

clear all; close all; clc;

%% User parameters
centerFreq = 96.0e6; %RMF FM
basebandFs = 160e3;
plutoGain  = 65;
frameLen   = 16384;
audioFs    = 48000;
doStereo   = true;          % set true to attempt stereo decode (L+R and L-R)
saveAudio  = true;           % set true to save .wav file
audioFilename = 'pluto_fm_output_plot.wav';
saveIQ     = false;          % save raw IQ to file (for later analysis)
iqFilename = 'pluto_fm_iq.dat';
deemphasisTau = 50e-6;       % de-emphasis time constant: Europe typically 50e-6, US 75e-6

%% Create Pluto receiver System object
rx = comm.SDRRxPluto( ...
    'CenterFrequency', centerFreq, ...
    'BasebandSampleRate', basebandFs, ...
    'OutputDataType','double', ...
    'RadioID','usb:0', ...          % adjust if you have multiple radios
    'GainSource','Manual', ...
    'Gain', plutoGain, ...
    'SamplesPerFrame', frameLen);

% Audio output
deviceWriter = audioDeviceWriter('SampleRate', audioFs);

% Spectrum Visualizer Setup
Nfft = 4096;
f = linspace(-basebandFs/2, basebandFs/2, Nfft)/1e6; % frequency axis in MHz
hFig = figure('Name','Live Spectrum','NumberTitle','off');
hAx = axes('Parent',hFig);
hPlot = plot(hAx, f, zeros(1,Nfft));
xlabel(hAx,'Frequency [MHz]');
ylabel(hAx,'Magnitude [dB]');
title(hAx,'Live Spectrum near tuned frequency');
grid(hAx,'on');
ylim(hAx, [-120 -20]);


% Prepare buffers / filters
% 1) FM demodulator: we'll use instantaneous phase differentiation
% 2) Lowpass for mono (15 kHz)
% 3) For stereo decode: bandpass for the 23-53 kHz region around stereo subcarrier

% Design filters (use dsp.LowpassFilter/dsp.BiquadFilter or designfilt)
monoLP = designfilt('lowpassiir', 'PassbandFrequency', 15e3, ...
    'StopbandFrequency', 18e3, 'SampleRate', basebandFs, 'PassbandRipple', 1, 'StopbandAttenuation', 60);

% For stereo path, need to isolate 23..53 kHz (L-R DSB around 38k) then shift 19k subcarrier
stereoBP = designfilt('bandpassiir', 'FilterOrder', 8, ...
    'HalfPowerFrequency1', 23e3, 'HalfPowerFrequency2', 53e3, 'SampleRate', basebandFs);

% De-emphasis: simple first-order IIR (digital) at audioFs (we will apply after resampling)
% Create continuous-time RC => discrete using bilinear transform later

% Output storage
audioBuf = [];
if saveIQ
    fidIQ = fopen(iqFilename,'wb');
end

disp('Starting FM reception. Press Ctrl-C to stop.');

% runtime variables for stereo demodulation
tframe = (0:frameLen-1)'/basebandFs; % local time vector for phase shift operations

% Let's set a fixed amount of time for recording
durationSec = 45; 
tStart = tic;

try
    while toc(tStart) < durationSec
        % Receive a frame of complex baseband IQ (Nx1 complex)
        rxIQ = rx();

        % --- Spectrum update ---
        X = fftshift(fft(rxIQ, Nfft));
        P = 20*log10(abs(X)/Nfft + 1e-12);
        if ishandle(hFig)
            set(hPlot, 'YData', P);
            drawnow limitrate;
        end


        if isempty(rxIQ)
            pause(0.01);
            continue;
        end

        % Optionally save IQ
        if saveIQ
            fwrite(fidIQ, [real(rxIQ) imag(rxIQ)] , 'double');
        end

        % FM demodulation: instantaneous frequency via angle(diff)
        % Compute angle difference: unwrap for stability
        ang = angle(rxIQ);
        % instantaneous frequency (rad/sample) -> convert to Hz by dividing by 2*pi and multiply by sample rate
        dAng = [0; diff(unwrap(ang))];                % rad
        instFreq = (basebandFs/(2*pi)) * dAng;       % Hz (relative instantaneous frequency)
        % The audio baseband (AF) is proportional to instFreq (remove DC)
        % Add this scaling after calculating instFreq:
        maxDeviation = 75000; % Standard maximum frequency deviation for WFM
        % Scale the audio baseband so that a 75 kHz deviation results in a level close to 1
        audioBase = instFreq / maxDeviation;

        % At this point audioBase contains baseband audio up to ~53 kHz (stereo includes subcarrier)
        % Filter to get mono (L+R)
        mono = filtfilt(monoLP, audioBase);  % zero-phase for clarity

        if doStereo
            % Isolate stereo subband (L-R around 38k DSB)
            stereoSub = filtfilt(stereoBP, audioBase);  % bandpassed around 23-53k
            % Multiply by 2*cos(2*pi*19k*t) to shift DSB suppressed carrier to baseband
            % For continuous stream we must maintain a phase that continues across frames.
            if ~exist('phaseOffset','var') || isempty(phaseOffset)
                phaseOffset = 0;
            end
            if isempty(phaseOffset), phaseOffset = 0; end
            % build time vector considering phaseOffset
            N = length(stereoSub);
            tvec = (0:N-1)'/basebandFs;
            carrier = 2*cos(2*pi*19e3*tvec + phaseOffset);
            demodLR = stereoSub .* carrier;  % now contains L-R around baseband
            % update phaseOffset for continuity
            phaseOffset = phaseOffset + 2*pi*19e3*(N/basebandFs);
            % lowpass to recover L-R (cutoff 15 kHz)
            LR = filtfilt(monoLP, demodLR);

            % compute left/right
            left  = 0.5*(mono + LR);
            right = 0.5*(mono - LR);

            % combine to stereo signal (two cols)
            % Now resample from basebandFs to audioFs (48k). Use resample to anti-alias.
            leftRes  = resample(left,  audioFs, basebandFs);
            rightRes = resample(right, audioFs, basebandFs);
            audioFrame = [leftRes, rightRes];

        else
            % Mono path: resample mono to audioFs
            monoRes = resample(mono, audioFs, basebandFs);
            % make stereo by duplicating channel (so playback uses 2 channels)
            audioFrame = [monoRes, monoRes];
        end

        % De-emphasis filter (first-order RC) applied at audioFs
        % y[n] = a * x[n] + (1-a) * y[n-1], where a = dt/(tau + dt), dt = 1/fs
        a = (1/audioFs) / (deemphasisTau + 1/audioFs);
        % apply across each channel with simple IIR
        audioFrame = filter(a, [1 (a-1)], audioFrame);  % single-pole lowpass style de-emphasis

        % Play audio
        deviceWriter(audioFrame);

        % Append to buffer for saving
        if saveAudio
            audioBuf = [audioBuf; audioFrame]; %#ok<AGROW>
        end
    end

catch ME
    % Clean up on Ctrl-C or error
    release(rx);
    release(deviceWriter);
    if saveIQ
        fclose(fidIQ);
    end
    warning('Reception stopped: %s', ME.message);
end

% Save audio
if saveAudio && ~isempty(audioBuf)
    % normalize and write to WAV
    audioBuf = audioBuf / max(abs(audioBuf(:)));
    audiowrite(audioFilename, audioBuf, audioFs);
    fprintf('Audio saved to %s\n', audioFilename);
end

disp('Done.');
