% Zadanie 2.
% Wielokrokowa metoda predyktor–korektor Adamsa 4. rzêdu ze sta³ym krokiem.


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


% Kroki calkowania (odpowiadajace kolejnym warunkom poczatkowym)
h = [0.004 0.004 0.005 0.2];


% Zbyt duze kroki calkowania (odpowiadajace kolejnym warunkom poczatkowym)
h_big = [0.005 0.008 0.01 0.3];


% Liczba warunkow poczatkowych
L = length(h);


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
    
    % Obliczenie rozwiazan i bledow oszacowan dla "dobrego" kroku h
    tic;
    [y, delta_h] = Adams4_const_h(f, m, a, b, h(i), y0(i,:));
    T(i) = toc;
    
    
    % Wykres fazowy dla "dobrego" kroku h
    figure(q);
    hold on;
    plot( y(:,1), y(:,2) );
    grid on;
    title( ['Rozwiazanie w przestrzeni fazowej - warunki poczatkowe nr ', num2str(i), ', krok h=', num2str( h(i) )] );
    xlabel('x1');
    ylabel('x2');
    hold off;
    q = q+1;
    
    
    % Obliczenie rozwiazan i bledow oszacowan (zbyt duzy krok h)
    [y_big, delta_h_big] = Adams4_const_h(f, m, a, b, h_big(i), y0(i,:));
    
    
    % Wykres fazowy dla zbyt "duzego" kroku h
    figure(q);
    hold on;
    plot( y_big(:,1), y_big(:,2) );
    grid on;
    title( ['Rozwiazanie w przestrzeni fazowej - warunki poczatkowe nr ', num2str(i), ', zbyt "duzy" krok h=', num2str( h_big(i) )] );
    xlabel('x1');
    ylabel('x2');
    hold off;
    q = q+1;
    
    
    % Wyznaczenie wektora (osi poziomej) do wykresow bledow aproksymacji
    x = a:h(i):b;
    
    
    % Wyznaczenie liczby iteracji
    iter(i) = length(x);
    
    
    % Wykres modulow bledow aproksymacji - zmienna zalezna x1 ("dobry" krok h)
    figure(q);
    hold on;
    plot( x([2:length(x)]), abs( delta_h(:,1) ) );
    grid on;
    title( ['Oszacowanie modulow bledow aproksymacji - zmienna zalezna x1 - warunki poczatkowe nr ', num2str(i), ', "dobry "krok h=', num2str( h(i) )] );
    xlabel('t');
    ylabel('Oszacowanie bledu aproksymacji');
    hold off;
    q = q+1;
    
    
    % Wykres bledu aproksymacji - zmienna zalezna x2 ("dobry" krok h)
    figure(q);
    hold on;
    plot( x([2:length(x)]), abs( delta_h(:,2) ) );
    grid on;
    title( ['Oszacowanie modulow bledow aproksymacji - zmienna zalezna x2 - warunki poczatkowe nr ', num2str(i), ', "dobry" krok h=', num2str( h(i) )] );
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
    plot( x([2:length(x)]), delta_h_aggr );
    grid on;
    title( ['Norma maksimum bledu aproksymacji - warunki poczatkowe nr ', num2str(i), ', "dobry" krok h=', num2str( h(i) )] );
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
    title( ['Porownanie rozwiazan - przestrzen fazowa - warunki poczatkowe nr ', num2str(i), ', "dobry" krok h=', num2str( h(i) )] );
    xlabel('x1');
    ylabel('x2');
    legend('Adams4 const h', 'ode45', 'Location', 'North');
    hold off;
    q = q+1;
    
    
    % Zmienna zalezna x1
    figure(q);
    hold on;
    plot( x, y(:,1), x_ode, y_ode(:,1) );
    grid on;
    title( ['Porownanie rozwiazan - zmienna zalezna x1 - warunki poczatkowe nr ', num2str(i), ', "dobry" krok h=', num2str( h(i) )] );
    xlabel('t');
    ylabel('x1');
    legend('Adams4 const h', 'ode45');
    hold off;
    q = q+1;
    
    
    % Zmienna zalezna x2
    figure(q);
    hold on;
    plot( x, y(:,2), x_ode, y_ode(:,2) );
    grid on;
    title( ['Porownanie rozwiazan - zmienna zalezna x2 - warunki poczatkowe nr ', num2str(i), ', "dobry" krok h=', num2str( h(i) )] );
    xlabel('t');
    ylabel('x2');
    legend('Adams4 const h', 'ode45');
    hold off;
    q = q+1;
    
end

