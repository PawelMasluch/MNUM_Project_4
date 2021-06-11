function [x, y, delta_h, error_code] = RK4_var_h(f, m, a, b, h, h_min, y0, eps_w, eps_b)
% Metoda Rungego-Kutty RK4 ze zmiennym krokiem,
% sluzaca do rozwiazywania rownan 
% i ukladow rownan rozniczkowych zwyczajnych.

% x - wektor wartosci zmiennej niezaleznej
% y - wartosci zmiennych zaleznych, 
%     dla okreslonego przedzialu wartosci zmiennej niezaleznej
% delta_h - oszacowania bledu pojedynczego kroku o d³ugoœci h
%           dla kazdej zmiennej zaleznej
% error_code - kod bledu (0 - brak bledow, 
%                         1 - niemo¿liwe rozwiazanie z zadana dokladnoscia h_min)
% f - funkcje prawych stron (w tablicy komorkowej)
% m - liczba rownan (i zmiennych zaleznych zarazem)
% [a, b] - przedzial zmiennej niezaleznej, 
%          dla ktorego poszukujemy rozwiazania
% y0 - wektor warunkow poczatkowych
% h - poczatkowy krok calkowania
% h_min - minimalny krok calkowania
% eps_w - dokladnosc wzgledna
% eps_b - dokladnosc bezwzgledna


% Poczatkowo zakladamy, ze nie bedzie bledow
error_code = 0;


% Wspolczynnik bezpieczenstwa
% (do wyznaczania nowego kroku)
s = 0.9;


% Heurystyczny wspolczynnik maksymalnego wzrostu kolejnego kroku
beta = 5;


% Rzad metody
p = 4;


% Wartosci pomocnicze k, charakterystyczne dla naszej metody;
% wartosci te wyliczamy w ramach konkretnej iteracji
k = zeros(4, m);


% Wartosci zmiennej niezaleznej, 
% dla ktorych obliczymy wartosci zmiennych zaleznych
x(1) = a;


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


% Licznik iteracji calkowania
nr = 2;


% Pomocnicza flaga - czy zakonczyc caly proces calkowania (0-nie, 1-tak)
koniec = 0;


% Dopoki iterujemy, dopoty nie osiagnelismy konca przedzialu 
% zmiennej niezaleznej
while koniec == 0
    
    % Wektor zmiennych zaleznych, obliczanych poprzez 2 iteracje
    % z krokiem h/2 
    y_tmp = y(nr-1,:);


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


    % Obliczenie oszacowania bledu pojedynczego kroku o dlugosci h
    delta_h(nr-1, :) = (2^p / (2^p - 1)) * ( y_tmp - y(nr,:) );


    % Wektor parametrow dokladnosci obliczen
    eps = abs( y(nr, :) ) * eps_w   +   eps_b;


    % Obliczenie wspolczynnika modyfikacji kroku
    alfa = (   eps ./ abs( delta_h(nr-1, :) )   ) .^ ( 1/(p+1) );
    alfa = min(alfa);


    % Proponowana korekta dlugosci kroku
    h_tmp = s*alfa*h;


    if s*alfa >= 1 % jesli h_tmp >= h
        if x(nr-1) + h > b
            y = y( [1:length(y)-1], : ); % pozbywamy sie ostatniego (nadmiarowego) wiersza
            koniec = 1; % koniec calkowania
        else
            x(nr) = x(nr-1) + h; % nowa wartosc zmiennej niezaleznej
            h = min( [h_tmp   beta*h   b - x(nr-1)] ); % kolejny krok
            nr = nr + 1; % kolejna iteracja
        end
    else
        if h_tmp < h_min
            x = [];
            y = [];
            delta_h = [];
            error_code = 1; % niemo¿liwe rozwiazanie z zadana dokladnoscia h_min
            koniec = 1; % koniec calkowania
        else
            h = h_tmp;
        end
    end
    
end

end


