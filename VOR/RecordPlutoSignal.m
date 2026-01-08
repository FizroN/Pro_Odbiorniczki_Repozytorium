function [receivedSignal, metadata] = RecordPlutoSignal(filename, centerFreq, sampleRate, duration, gain)
%RECORDPLUTOSIGNAL Nagrywa sygnał z PlutoSDR z detekcją utraconych próbek.
%   Zapisuje wynik do pliku .mat i zwraca dane do workspace.

    % 1. Konfiguracja domyślna
    if nargin < 5
        gain = []; 
    end
    
    fprintf('Inicjalizacja PlutoSDR (Fc: %.2f MHz, Fs: %.2f MHz)...\n', centerFreq/1e6, sampleRate/1e6);
    
    try
        % Konfiguracja obiektu odbiornika
        rx = sdrrx('Pluto');
        rx.CenterFrequency = centerFreq;
        rx.BasebandSampleRate = sampleRate;
        rx.SamplesPerFrame = 65536; % Duży bufor dla stabilności USB
        rx.OutputDataType = 'double';
        
        if isempty(gain)
            rx.GainSource = 'AGC Fast Attack';
        else
            rx.GainSource = 'Manual';
            rx.Gain = gain;
        end
        
    catch ME
        error('Błąd połączenia z PlutoSDR: %s', ME.message);
    end
    
    % 2. Obliczenia parametrów nagrywania
    samplesPerFrame = rx.SamplesPerFrame;
    totalSamples = floor(duration * sampleRate);
    numFrames = ceil(totalSamples / samplesPerFrame);
    
    % Prealokacja pamięci (Szybkość!)
    receivedSignal = complex(zeros(numFrames * samplesPerFrame, 1));
    
    % 3. Rozgrzewka (Wyrzucenie stanów nieustalonych AGC/DC)
    % Pobieramy kilka ramek i wyrzucamy do kosza
    for k = 1:20 
        rx();
    end
    
    fprintf('Nagrywanie sygnału... ');
    
    % 4. Pętla Nagrywania z detekcją błędów
    lostSamplesCount = 0;
    tic;
    for k = 1:numFrames
        % Pobranie ramki danych
        [dataFrame, len, lost] = rx();
        
        % Detekcja zerwania strumienia
        if lost > 0
            lostSamplesCount = lostSamplesCount + 1;
            fprintf('!'); % Wykrzyknik w konsoli oznacza zgubione próbki
        end
        
        % Wklejenie danych do głównej tablicy
        startIndex = (k-1) * samplesPerFrame + 1;
        endIndex = k * samplesPerFrame;
        receivedSignal(startIndex:endIndex) = dataFrame;
    end
    actualDuration = toc;
    fprintf('\n');
    
    % 5. Raportowanie błędów
    if lostSamplesCount > 0
        warning('UWAGA: Zgubiono próbki w %d ramkach! Sygnał nieciągły.', lostSamplesCount);
        fprintf('Sugerowane: Zmniejsz SampleRate lub zamknij inne aplikacje.\n');
    else
        fprintf('Sukces: Brak zgubionych próbek (Sygnał ciągły).\n');
    end
    
    % Przycięcie nadmiaru bufora i zwolnienie radia
    receivedSignal = receivedSignal(1:totalSamples);
    release(rx);
    
    % 6. Zapis metadanych i pliku
    metadata.centerFreq = centerFreq;
    metadata.sampleRate = sampleRate;
    metadata.duration = duration;
    metadata.timestamp = datetime('now');
    metadata.gainUsed = rx.Gain;
    metadata.lostSamples = lostSamplesCount;
    
    % Zapisuj tylko wtedy, gdy filename nie jest pusty
    if ~isempty(filename)
        save(filename, 'receivedSignal', 'metadata');
        fprintf('Zapisano do pliku .mat: %s\n', filename);
    else
        fprintf('Pominięto zapis pliku .mat (zgodnie z ustawieniem).\n');
    end
    
end