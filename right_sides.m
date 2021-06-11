function dy = right_sides(x,y)
% Funkcja obliczajaca funkcje prawych stron 
% dla ukladu rownan rozniczkowych z zadania
% (na potrzeby metody ode45).

% x - zmienna niezalezna
% y - wektor zmiennych zaleznych

dy = zeros(2,1);
dy(1) =  y(2) + y(1)*(0.9 - y(1)^2 - y(2)^2);
dy(2) = -y(1) + y(2)*(0.9 - y(1)^2 - y(2)^2);

end

