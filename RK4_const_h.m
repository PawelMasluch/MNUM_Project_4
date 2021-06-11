function [y, delta_h] = RK4_const_h(f, m, a, b, h, y0)
% Metoda Rungego-Kutty RK4 ze stalym krokiem,
% sluzaca do rozwiazywania rownan 
% i ukladow rownan rozniczkowych zwyczajnych.

% y - wartosci zmiennych zaleznych, 
%     dla okreslonego przedzialu wartosci zmiennej niezaleznej
% delta_h - oszacowania bledu pojedynczego kroku o d³ugoœci h
%           dla kazdej zmiennej zaleznej
% f - funkcje prawych stron (w tablicy komorkowej)
% m - liczba rownan (i zmiennych zaleznych zarazem)
% [a, b] - przedzial zmiennej niezaleznej, 
%          dla ktorego poszukujemy rozwiazania
% y0 - wektor warunkow poczatkowych
% h - krok calkowania


% Rzad metody
p = 4;


% Wartosci zmiennej niezaleznej, 
% dla ktorych obliczymy wartosci zmiennych zaleznych
x = a:h:b;


% Liczba krokow calkowania
n = length(x);


% Wartosci pomocnicze k, charakterystyczne dla naszej metody;
% wartosci te wyliczamy w ramach konkretnej iteracji
k = zeros(4, m);


% y - wartosci zmiennych zaleznych:
% wiersze - wartosci wszystkich zmiennych zaleznych 
%           dla danej iteracji, 
% kolumny - wartosci okreslonych zmiennych zaleznych 
%           we wszystkich iteracjach;
% Inicjacja warunkami poczatkowymi
y(1, :) = y0;


% Oszacowania bledu pojedynczego kroku o d³ugoœci h
% dla kazdej zmiennej zaleznej
delta_h(1,:) = zeros(m, 1);


% (Prawie) kazda iteracja
for nr=2:n
    
    % Wektor zmiennych zaleznych, obliczanych poprzez 2 iteracje
    % z krokiem h/2 
    y_tmp = y(nr-1,:);
    
    % ----------------------------------
    
    % 2 iteracje z krokiem h/2
    % (celem oszacowania b³êdu pojedynczego kroku o d³ugoœci h)
    for j=1:2
        
        % Wartosci k1 dla kazdej zmiennej zaleznej
        for i=1:m
            k(1,i) = feval( f{i}, x(nr-1), y_tmp );
        end


        % Wartosci k2 dla kazdej zmiennej zaleznej
        for i=1:m
            k(2,i) = feval( f{i}, x(nr-1) + j*h/4,   y_tmp + j*h/4 * k(1,:) );
        end


        % Wartosci k3 dla kazdej zmiennej zaleznej
        for i=1:m
            k(3,i) = feval( f{i},   x(nr-1) + j*h/4,   y_tmp + j*h/4 * k(2,:) );
        end


        % Wartosci k4 dla kazdej zmiennej zaleznej
        for i=1:m
            k(4,i) = feval( f{i},   x(nr-1) + j*h/2,   y_tmp + j*h/2 * k(3,:) );
        end


        % Obliczenie wartosci zmiennych zaleznych w aktualnym kroku
        y_tmp = y_tmp + h/12 * ( k(1,:) + 2*k(2,:) + 2*k(3,:) + k(4,:) );
    end
    
    % --------------------------------------
    
    
    % 1 iteracja z krokiem h
    
    % Wartosci k1 dla kazdej zmiennej zaleznej
    for i=1:m
        k(1,i) = feval( f{i}, x(nr-1), y(nr-1,:) );
    end
    
    
    % Wartosci k2 dla kazdej zmiennej zaleznej
    for i=1:m
        k(2,i) = feval( f{i}, x(nr-1) + h/2,   y(nr-1,:) + h/2 * k(1,:) );
    end
    
    
    % Wartosci k3 dla kazdej zmiennej zaleznej
    for i=1:m
        k(3,i) = feval( f{i},   x(nr-1) + h/2,   y(nr-1,:) + h/2 * k(2,:) );
    end
    
    
    % Wartosci k4 dla kazdej zmiennej zaleznej
    for i=1:m
        k(4,i) = feval( f{i},   x(nr-1) + h,   y(nr-1,:) + h * k(3,:) );
    end
    
    
    % Obliczenie wartosci zmiennych zaleznych w aktualnym kroku
    y(nr,:) = y(nr-1,:) + h/6 * ( k(1,:) + 2*k(2,:) + 2*k(3,:) + k(4,:) );
    
    % --------------------------
    
    
    % Obliczenie oszacowania bledu pojedynczego kroku o dlugosci h
    delta_h(nr-1, :) = (2^p / (2^p - 1)) * ( y_tmp - y(nr,:) );
    
end


end

