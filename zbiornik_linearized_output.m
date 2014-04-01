function dydt = zbiornik_linearized_output(t, x)
% zbiornik_linearized_output( t, x, u ) - funkcja obliczajaca wyjscie
% zlinearyzowanego modelu w postaci rownan stanu dla ukladu dynamicznego
% zbiornika z mieszaniem.
%
%   ARGUMENTY:
%     t - czas
%     x - wartosc stanu
%   WARTOSCI WYJSCIOWE:
%     dydt - wartosc pochodnej po czasie funckji stanu
% 
% Funkcja korzysta ponadto ze zmiennych globalnych:
% 
%     C
% 
% ktora jest macierza opisujaca zaleznosc wyjscia od stanu obiektu. Musi
% ona zostac zainicjalizowana przed wywolaniem funkcji.

global C;

dydt = [ ...
   C(1,1)*x(1) + C(1,2)*x(2);
   C(2,1)*x(1) + C(2,2)*x(2);
];
