% Fragment kodu dla Zadania 21.6 - Modyfikacja NLMS na LMS
% Zakładamy, że zmienne u, v, M (długość filtra) są już zdefiniowane jak w Listingu 21.5

% Inicjalizacja
bv = zeros(1, M);      % Bufor wejściowy
g = zeros(1, M);       % Wagi korektora (początkowo zera)
ghist = [];            % Do wykresu historii wag

% WAŻNE: Dla zwykłego LMS współczynnik 'mi' musi być znacznie mniejszy niż dla NLMS!
% W NLMS było np. 0.1, dla LMS przy silnym sygnale spróbuj rzędu 0.001 - 0.01
mi_LMS = 0.005; 

% Pętla adaptacyjna
for n = 1:length(u)
    % 1. Aktualizacja bufora wejściowego (przesunięcie + nowa próbka)
    % v to sygnał odebrany (zniekształcony)
    bv = [v(n) bv(1:M-1)]; 
    
    % 2. Obliczenie wyjścia filtru (estymata)
    uest(n) = sum(bv .* conj(g)); % splot wejścia z wagami
    
    % 3. Obliczenie błędu
    % u(n) to sygnał referencyjny (pilot/nagłówek)
    err(n) = u(n) - uest(n);
    
    % 4. Aktualizacja wag (Tu jest zmiana względem NLMS!)
    
    % Wersja NLMS (była w książce):
    % energy = (bv * bv'); 
    % g = g + mi * conj(err(n)) * bv / (energy + 1e-10);
    
    % Wersja LMS (ZADANIE):
    % Brak normalizacji przez energię.
    g = g + mi_LMS * conj(err(n)) * bv;
    
    % Zapis historii pierwszej wagi do wykresu
    ghist = [ghist g(1)];
end

% Wyświetlenie wyników (jak w oryginale)
figure; plot(abs(ghist)); title('Adaptacja wagi g(1) (LMS)'); grid on;
figure; 
plot(real(uest(end-100:end)), 'r'); hold on; 
plot(real(u(end-100:end)), 'b--'); 
title('Porównanie: Czerwony=Wyjście korektora, Niebieski=Oczekiwane');
legend('Equalized', 'Desired');
grid on;