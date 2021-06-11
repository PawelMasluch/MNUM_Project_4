function [ y, delta ] = Adams4_const_h( f, m, a, b, h, y0 )
% Wielokrokowa metoda predyktor–korektor Adamsa 4. rzêdu ze sta³ym krokiem,
% sluzaca do rozwiazywania rownan 
% i ukladow rownan rozniczkowych zwyczajnych.
% Metoda predyktora: metoda (jawna)   4-krokowa Adamsa-Bashfortha (4. rzedu)
% Metoda korektora: metoda (niejawna) 3-krokowa Adamsa-Moultona (4. rzedu)
% Schemat metody: P4 E K3 E  

% y - wartosci zmiennych zaleznych, 
%     dla okreslonego przedzialu wartosci zmiennej niezaleznej
% delta - oszacowania bledu pojedynczego kroku o d³ugoœci h
%           dla kazdej zmiennej zaleznej
% f - funkcje prawych stron (w tablicy komorkowej)
% m - liczba rownan (i zmiennych zaleznych zarazem)
% [a, b] - przedzial zmiennej niezaleznej, 
%          dla ktorego poszukujemy rozwiazania
% y0 - wektor warunkow poczatkowych
% h - krok calkowania



% Liczba krokow predyktora
k_pred = 4;


% Liczba krokow korektora
k_corr = 3;


% Wspolczynniki dla metody predyktora 
% (tj. metody (jawnej) 4-krokowej Adamsa-Bashfortha)
beta_pred = 1/24 * [ 55 -59 37 -9 ];


% Wspolczynniki dla metody korektora
% (tj. metody (niejawnej) 3-krokowej Adamsa-Moultona)
beta_corr = 1/24 * [ 9 19 -5 1 ];


% Wartosci zmiennej niezaleznej, 
% dla ktorych obliczymy wartosci zmiennych zaleznych
x = a:h:b;


% Liczba krokow calkowania
n = length(x);


% y - wartosci zmiennych zaleznych:
% wiersze - wartosci wszystkich zmiennych zaleznych 
%           dla danej iteracji, 
% kolumny - wartosci okreslonych zmiennych zaleznych 
%           we wszystkich iteracjach;
% Inicjacja warunkami poczatkowymi
y(1, :) = y0;


% Oszacowania bledu pojedynczego kroku o dlugosci h
% dla kazdej zmiennej zaleznej
delta(1,:) = zeros(m, 1);


% (Prawie) kazda iteracja
for nr=2:n
    
    % 1) Predykcja
    
    % Wektor wartosci predykcji
    y_pred = zeros(1,m);
    
    % Kazda zmienna zalezna
    for i=1:m
        
        % Obliczenie sumy pomocniczej 
        sum = 0;
        for j=1:k_pred
            if nr - j >= 1
                sum = sum + beta_pred(j) * feval( f{i}, x(nr-j), y(nr-j, :) );
            end
        end
        
        % Obliczenie wartosci predykcji i-tej zmiennej zaleznej
        % w aktualnej iteracji
        y_pred(i) = y(nr-1, i) + h * sum;
    end
    
    
    % 2) Ewaluacja
    
    % Wektor wartosci funkcji prawych stron
    % dla wartosci predykcji
    f_eval = zeros(1,m);
    
    for i=1:m
        f_eval(i) = feval( f{i}, x(nr), y_pred );
    end
    
    
    % 3) Korekcja
    
    % Kazda zmienna zalezna
    for i=1:m
        
        % Obliczenie sumy pomocniczej 
        sum = 0;
        for j=2:k_corr+1
            if nr - j >= 1
                sum = sum + beta_corr(j) * feval( f{i}, x(nr-j+1), y(nr-j+1, :) );
            end
        end
        
        % Obliczenie wartosci i-tej zmiennej zaleznej
        % w aktualnej iteracji
        y(nr, i)   =   y(nr-1, i)   +   h*sum   +   h * beta_corr(1) * f_eval(i);
    end
    
    
    % 4) Ewaluacja - nie jest konieczna
    
    
    % Obliczenie wektora estymat bledow aproksymacji
    delta(nr-1, :) = -19/270 * ( y_pred - y(nr, :) );
    
end



end

