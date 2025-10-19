% ==========================================================
% SDR_IQ_Recorder.m — Real-Time Spectrum View using ADALM-PLUTO
% ==========================================================
clear; close all; clc;

%% --- RECORDING PARAMETERS ---
fc = 2422;                % carrier frequency in MHz
BW = 25;               % bandwidth in MHz
SamplesPerFrame = 76800;  % number of IQ samples per frame

%% --- SDR DEVICE INITIALIZATION ---
sdrdev                         % show installed SDR devices
rxPluto = sdrrx('Pluto', ...
    'RadioID','usb:0', ...
    'CenterFrequency', fc * 1e6, ...
    'BasebandSampleRate', BW * 1e6, ...
    'SamplesPerFrame', SamplesPerFrame, ...
    'OutputDataType','double');

%% --- REAL-TIME SPECTRUM VIEW SETUP ---
fs = BW * 1e6;
Nfft = 4096;
f = linspace(-fs/2, fs/2, Nfft)/1e6; % frequency axis [MHz]

hFig = figure('Name','Live Spectrum','NumberTitle','off');
hAx = axes('Parent',hFig);
hPlot = plot(hAx, f, zeros(1,Nfft));
xlabel(hAx,'Frequency [MHz]');
ylabel(hAx,'Magnitude [dB]');
title(hAx,'Real-Time Spectrum from ADALM-PLUTO');
grid(hAx,'on');
ylim(hAx, [-120 -20]);

%% --- CONTINUOUS RECEIVE AND DISPLAY LOOP ---
disp('▶ Receiving and displaying live spectrum... Press CTRL+C to stop.');

while ishandle(hFig)
    [data, valid, overflow] = rxPluto();
    if valid
        X = fftshift(fft(data, Nfft));
        P = 20*log10(abs(X)/Nfft + 1e-12);  % avoid log(0)
        set(hPlot, 'YData', P);
        drawnow limitrate;
    end
    if overflow
        warning('⚠ Overflow detected — samples dropped.');
    end
end
