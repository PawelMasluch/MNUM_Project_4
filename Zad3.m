% Zadanie 3.
% Metoda Rungego-Kutty RK4 ze zmiennym krokiem.


clear;
clc;


% Liczba rownan rozniczkowych (i zmiennych zaleznych zarazem)
m = 2;


% f - tablica komorkowa, zawierajaca uchwyty 
% do prawych stron rownan rozniczkowych
f = cell(m,1);


% Wstawienie (do tablicy f) uchwytow 
% do prawych stron rownan rozniczkowych
f{1} = @(x, y)( y(2) + y(1)*(0.9 - y(1)^2 - y(2)^2) );
f{2} = @(x, y)( -y(1) + y(2)*(0.9 - y(1)^2 - y(2)^2) );


% Przedzial [a,b] zmiennej niezaleznej, 
% dla ktorego bedziemy calkowali 
% uklad rownan rozniczkowych
a = 0;
b = 20;


% Warunki poczatkowe (kolejne wiersze)
y0 = [10    8;
       0    9; 
       8    0; 
    1e-3 1e-3];


% Poczatkoe kroki calkowania (odpowiadajace kolejnym warunkom poczatkowym)
h0 = [0.01465 0.025 0.02 0.2];


% Minimalne kroki calkowania
h_min = [1e-5 5e-5 5e-5 5e-5];


% Dokladnosc wzgledna
eps_w = 1e-5;


% Dokladnosc bezwzgledna
eps_b = 1e-5;


% Liczba warunkow poczatkowych
L = 4;


% Numery zestawow warunkow poczatkowych
nr = 1:L;


% Czasy dzialania solwerow
T = zeros(L,1);


% Liczby iteracji
iter = zeros(L,1);


% Licznik wykresow
q = 1;


% Wyznaczenie rozwiazan 
% dla roznych warunkow poczatkowych
for i=1:L
    
    % Obliczenie rozwiazan i bledow oszacowan
    tic;
    [x, y, delta_h, error_code] = RK4_var_h(f, m, a, b, h0(i), h_min(i), y0(i,:), eps_w, eps_b);
    T(i) = toc;
    
    
    % Wyznaczenie liczby iteracji
    iter(i) = length(x);
    
    
    if error_code == 0 % pomyslne wykonanie metody
        
        % Wykres modulow bledow aproksymacji - zmienna zalezna x1 ("dobry" krok h)
        figure(q);
        hold on;
        plot( x, abs( delta_h(:,1) ) );
        grid on;
        title( ['Oszacowanie modulow bledow aproksymacji - zmienna zalezna x1 - warunki poczatkowe nr ', num2str(i)] );
        xlabel('t');
        ylabel('Oszacowanie bledu aproksymacji');
        hold off;
        q = q+1;


        % Wykres bledu aproksymacji - zmienna zalezna x2 ("dobry" krok h)
        figure(q);
        hold on;
        plot( x, abs( delta_h(:,2) ) );
        grid on;
        title( ['Oszacowanie modulow bledow aproksymacji - zmienna zalezna x2 - warunki poczatkowe nr ', num2str(i)] );
        xlabel('t');
        ylabel('Oszacowanie bledu aproksymacji');
        hold off;
        q = q+1;


        % Wykres normy maksimum bledu aproksymacji ("dobry" krok h)
        delta_h_aggr = zeros( length(delta_h), 1 );
        for c=1:length(delta_h)
            delta_h_aggr(c) = norm( delta_h(c, :), Inf );
        end

        figure(q);
        hold on;
        plot( x, delta_h_aggr );
        grid on;
        title( ['Norma maksimum bledu aproksymacji - warunki poczatkowe nr ', num2str(i)] );
        xlabel('t');
        ylabel('Oszacowanie bledu aproksymacji');
        hold off;
        q = q+1;


        % Rozwiazanie przy pomocy polecenia ode45
        [x_ode, y_ode] = ode45(@right_sides,[a b], y0(i,:));


        % Wykresy porownawcze

        % Wykres fazowy
        figure(q);
        hold on;
        plot( y(:,1), y(:,2), y_ode(:,1), y_ode(:,2) );
        grid on;
        title( ['Porownanie rozwiazan - przestrzen fazowa - warunki poczatkowe nr ', num2str(i)] );
        xlabel('x1');
        ylabel('x2');
        legend('RK4 var h', 'ode45', 'Location', 'North');
        hold off;
        q = q+1;


        % Zmienna zalezna x1
        figure(q);
        hold on;
        plot( x, y(:,1), x_ode, y_ode(:,1) );
        grid on;
        title( ['Porownanie rozwiazan - zmienna zalezna x1 - warunki poczatkowe nr ', num2str(i)] );
        xlabel('t');
        ylabel('x1');
        legend('RK4 var h', 'ode45');
        hold off;
        q = q+1;


        % Zmienna zalezna x2
        figure(q);
        hold on;
        plot( x, y(:,2), x_ode, y_ode(:,2) );
        grid on;
        title( ['Porownanie rozwiazan - zmienna zalezna x2 - warunki poczatkowe nr ', num2str(i)] );
        xlabel('t');
        ylabel('x2');
        legend('RK4 var h', 'ode45');
        hold off;
        q = q+1;
         
    else % metoda nie zakonczyla sie pomyslnie
        
        fprintf( 'Dla warunkoow poczatkowych nr %s, metoda zakonczyla sie niepowodzeniem.\n', num2str(i) );
    end
    
end

