function dxdt = zbiornik(t, x, u, z)
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

global plant_C;
global plant_alpha;
global plant_T_C0;
global plant_T_H0;

dxdt = [ ...
   u(1) + u(2) + z(1) - plant_alpha*(x(1)/plant_C)^(1/6);
   plant_T_H0*u(1)/x(1) + plant_T_C0*u(2)/x(1) + z(1)*z(2)/x(1) - (x(2)/x(1))*(u(1) + u(2) + z(1))
];
