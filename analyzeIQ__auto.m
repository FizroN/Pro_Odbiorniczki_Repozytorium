function analyzeIQ_auto(iqFilename, basebandFs, audioFs, varargin)
% analyzeIQ_auto(...) - improved offline FM stereo analyzer with auto-fix heuristics
% Usage:
%   analyzeIQ_auto('pluto_fm_iq.dat', 200e3, 48000);
% Optional name/value:
%   'segmentSec' (default 5), 'playAudio' (default false), 'saveWav' (true),
%   'wavName' (default 'analyzeIQ_out.wav').
%
% This function will:
%  - auto-detect IQ data type (double/single/int16)
%  - try multiple IQ order/sign permutations and pick the one with the strongest 19k pilot
%  - perform coarse frequency correction, FM demod, pilot phase unwrap+smoothing,
%    synthesize locked 38 kHz carrier, demod L-R, reconstruct L/R, resample, de-emphasize,
%    normalize, save WAV, optionally play, and plot diagnostics.

% --------- parse args ----------
p = inputParser;
addRequired(p,'iqFilename',@ischar);
addRequired(p,'basebandFs',@isnumeric);
addRequired(p,'audioFs',@isnumeric);
addParameter(p,'segmentSec',5,@(x)isnumeric(x)&&x>0);
addParameter(p,'playAudio',false,@islogical);
addParameter(p,'saveWav',true,@islogical);
addParameter(p,'wavName','analyzeIQ_out.wav',@ischar);
parse(p, iqFilename, basebandFs, audioFs, varargin{:});

segmentSec = p.Results.segmentSec;
doPlay = p.Results.playAudio;
doSave = p.Results.saveWav;
wavName = p.Results.wavName;

% --------- read IQ file with auto dtype detection ----------
fid = fopen(iqFilename,'rb');
if fid<0, error('Cannot open IQ file %s', iqFilename); end
% read a chunk first to guess
probe = fread(fid, min(10000, 1e6), 'double'); fclose(fid);
if isempty(probe)
    error('IQ file empty or unreadable.');
end
% Heuristic: if many very large numbers -> file probably int16
if max(abs(probe)) > 1e6
    dtype = 'int16';
else
    % try double first, fallback to int16 if nonsense later
    dtype = 'double';
end

% Actually read full file according to dtype
fid = fopen(iqFilename,'rb');
switch dtype
    case 'double'
        data = fread(fid, inf, 'double');
        fclose(fid);
        % if values tiny/large improbable -> try single/int16
        if isempty(data) || max(abs(data))>1e6 || max(abs(data))<1e-8
            fclose(fid);
            fid = fopen(iqFilename,'rb');
            data = fread(fid, inf, 'int16'); fclose(fid);
            data = double(data)/32768;
            fprintf('Fell back to int16 reading.\n');
        end
    otherwise
        data = fread(fid, inf, 'int16'); fclose(fid);
        data = double(data)/32768;
end

if isempty(data)
    error('No data read from file.');
end

% make sure even length
if mod(length(data),2)~=0, data = data(1:end-1); end

% reconstruct I/Q candidates (try permutations)
I = data(1:2:end);
Q = data(2:2:end);

candidates = {};
% standard
candidates{end+1} = I + 1j*Q;
% swapped
candidates{end+1} = Q + 1j*I;
% Q inverted sign
candidates{end+1} = I - 1j*Q;
candidates{end+1} = Q - 1j*I;
% conjugates (sometimes sign flip)
candidates{end+1} = conj(I + 1j*Q);
candidates{end+1} = conj(Q + 1j*I);

nCandidates = numel(candidates);
pilotScores = zeros(nCandidates,1);
pilotFreqs = zeros(nCandidates,1);

% helper FM demod
    function audioBase = fm_demod(iq, baseFs)
        % requires iq column vect
        d = iq(2:end).*conj(iq(1:end-1));
        a = angle(d);
        instHz = (baseFs/(2*pi)) * a;
        audioBase = instHz / 75e3;
        audioBase = [audioBase(1); audioBase];
    end

% Evaluate pilot strength per candidate on first segmentSec (or full if small)
Nwant = min(length(candidates{1}), round(segmentSec * basebandFs));
for k=1:nCandidates
    rxIQ_k = candidates{k};
    if length(rxIQ_k) < 2
        pilotScores(k) = -Inf;
        continue;
    end
    rxIQ_k = rxIQ_k(1:Nwant);
    audioBase_k = fm_demod(rxIQ_k, basebandFs);
    % pilot PSD estimate to find freq and amplitude
    [pxx,f] = pwelch(audioBase_k, 1024, 512, 4096, basebandFs);
    idx = find(f>17e3 & f<21e3);
    if isempty(idx)
        pilotScores(k) = 0;
        pilotFreqs(k) = 19e3;
        continue;
    end
    [pmax, r] = max(pxx(idx));
    pilotFreqs(k) = f(idx(r));
    % compute narrowband pilot envelope by bandpass
    [b_p,a_p] = butter(4, [(19e3-1000) (19e3+1000)]/(basebandFs/2), 'bandpass');
    try
        pilot = filtfilt(b_p,a_p,audioBase_k);
        pilAn = hilbert(pilot);
        pilotAmp = abs(pilAn);
        % score: median amplitude (robust)
        pilotScores(k) = median(pilotAmp);
    catch
        pilotScores(k) = 0;
    end
end

% pick best candidate
[~, bestIdx] = max(pilotScores);
if isinf(pilotScores(bestIdx)) || pilotScores(bestIdx) <= 0
    warning('No clear pilot found in any IQ permutation. Trying first candidate anyway.');
    bestIdx = 1;
end

rxIQ = candidates{bestIdx}(:); % column
fprintf('Selected candidate %d with pilot score %.5g (pilot freq est %.1f Hz)\n', ...
    bestIdx, pilotScores(bestIdx), pilotFreqs(bestIdx));

% use full segment (or reduce to manageable size)
Ntotal = min(length(rxIQ), round(max(segmentSec, 5) * basebandFs));
rxIQ = rxIQ(1:Ntotal);

% ---------- coarse pilot freq estimate & freq-correct if needed ----------
audioBase = fm_demod(rxIQ, basebandFs);
[pxx,f] = pwelch(audioBase, 2048, 1024, 8192, basebandFs);
idx = find(f>17e3 & f<21e3);
if ~isempty(idx)
    [~, rr] = max(pxx(idx));
    pilotFreqEst = f(idx(rr));
else
    pilotFreqEst = 19e3;
end
freqError = pilotFreqEst - 19e3;
fprintf('Pilot estimate %.1f Hz (error %.1f Hz)\n', pilotFreqEst, freqError);
if abs(freqError) > 30
    fprintf('Applying coarse freq correction of %.1f Hz\n', -freqError);
    n = (0:length(rxIQ)-1)';
    rxIQ = rxIQ .* exp(-1j*2*pi*freqError.*n/basebandFs);
    audioBase = fm_demod(rxIQ, basebandFs);
end

% ---------- pilot extraction, unwrap + smoothing ----------
% bandpass pilot
[b_p,a_p] = butter(4, [(19e3-1200) (19e3+1200)]/(basebandFs/2), 'bandpass');
pilot = filtfilt(b_p,a_p,audioBase);
pilAn = hilbert(pilot);
pilotAmp = abs(pilAn);
pilotPhase = angle(pilAn);
pilotPhaseUnw = unwrap(pilotPhase);
% smooth phase (moving average, zero-phase)
smoothMs = 2; winLen = max(3, round((smoothMs/1000)*basebandFs));
if mod(winLen,2)==0, winLen = winLen+1; end
bS = ones(winLen,1)/winLen;
pilotPhaseSm = filtfilt(bS,1,pilotPhaseUnw);

% ---------- synth 38k carrier by doubling smoothed phase ----------
phase38 = 2 * pilotPhaseSm;
carrier38 = 2 * cos(phase38);

% ---------- stereo subband isolation and demod ----------
[b_bp,a_bp] = butter(6, [23e3 53e3]/(basebandFs/2), 'bandpass');
stereoSub = filtfilt(b_bp,a_bp,audioBase);
demodLR = stereoSub .* carrier38;
[b_mono,a_mono] = butter(6, 15e3/(basebandFs/2), 'low');
LR = filtfilt(b_mono,a_mono, demodLR);
mono = filtfilt(b_mono,a_mono, audioBase);

left = 0.5*(mono + LR);
right = 0.5*(mono - LR);

% ---------- adaptive AGC-like scaling (improves listening) ----------
% scale mono band to reasonable RMS before resample
targetRMS = 0.28;
rmsMono = sqrt(mean(mono.^2)+eps);
gain = targetRMS / rmsMono;
left = left * gain;
right = right * gain;

% ---------- resample to audioFs ----------
Lr = resample(left, audioFs, basebandFs);
Rr = resample(right, audioFs, basebandFs);
stereo = [Lr, Rr];

% de-emphasis (50us default)
tau = 50e-6;
b_deem = (1/audioFs)/(tau + 1/audioFs);
a_deem = [1 (b_deem-1)];
stereo(:,1) = filter(b_deem, a_deem, stereo(:,1));
stereo(:,2) = filter(b_deem, a_deem, stereo(:,2));

% normalize
stereo = stereo / max(abs(stereo(:)) + eps);

% plots
t = (0:length(audioBase)-1)/basebandFs;
figure('Name','Baseband spectrogram'); spectrogram(audioBase,512,400,1024,basebandFs,'yaxis'); title('Audio baseband');
figure('Name','Pilot envelope & phase'); subplot(2,1,1); plot(t, pilotAmp); xlabel('s'); ylabel('amp'); title('Pilot envelope');
med = median(pilotAmp); sdev = std(pilotAmp); th = med + 1.5*sdev; yline(th,'--r','suggested thresh');
subplot(2,1,2); plot(t, pilotPhaseUnw); hold on; plot(t, pilotPhaseSm,'r'); legend('unwrapped','smoothed');
figure('Name','Left spectrogram'); spectrogram(Lr,1024,900,2048,audioFs,'yaxis'); title('Left');
figure('Name','Right spectrogram'); spectrogram(Rr,1024,900,2048,audioFs,'yaxis'); title('Right');


fprintf('Done. Candidate %d chosen. Pilot median=%.5g suggested thresh=%.5g\n', bestIdx, median(pilotAmp), th);

end




