function dydt = zbiornik_linearized_output(t, x)
% zbiornik( t, x, u ) - funkcja obliczajaca rownania stanu dla układu 
% dynamicznego zbiornika z mieszaniem
%
%   ARGUMENTY:
%     t - czas
%     x - wartosc stanu
%     u - wartosc sterowania
%   WARTOSCI WYJSCIOWE:
%     dxdt - wartosc pochodnej po czasie funckji stanu
% 
% Funkcja korzysta ponadto ze zmiennych globalnych:
% 
%     C, alpha, T_C0, T_H0, T_D0
% 
% które muszą zostać zainicjalizowane przed wywołaniem funkcji.

global A;
global B;
global C;

dydt = [ ...
   C(1,1)*x(1) + C(1,2)*x(2);
   C(2,1)*x(1) + C(2,2)*x(2);
];
