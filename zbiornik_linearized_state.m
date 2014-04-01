function dxdt = zbiornik_linearized_state(t, x, u)
% zbiornik_linearized_state( t, x, u ) - funkcja obliczjaca zlinearyzowane
% rownania stanu dla układu dynamicznego zbiornika z mieszaniem
%
%   ARGUMENTY:
%     t - czas
%     x - wartosc stanu
%     u - wejscia ( u1, u2 - sterowania; u3, u4 - zaklocenia)
%   WARTOSCI WYJSCIOWE:
%     dxdt - wartosc pochodnej po czasie funckji stanu
% 
% Funkcja korzysta ponadto ze zmiennych globalnych:
% 
%     A, B
% 
% bedacych macierzami opisujacymi zlinearyzowany model. Musza one zostac
% zainicjalizowane przed wywolaniem funkcji.

global A;
global B;

dxdt = [ ...
   A(1,1)*x(1) + A(1,2)*x(2) + B(1,1)*u(1) + B(1,2)*u(2) + B(1,3)*u(3) + B(1,4)*u(4);
   A(2,1)*x(1) + A(2,2)*x(2) + B(2,1)*u(1) + B(2,2)*u(2) + B(2,3)*u(3) + B(2,4)*u(4);
];
