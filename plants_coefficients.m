% Definicja stałych. 
% Są zadeklarowane jako zmienne globalne na rzecz wywoływania przez funkcję
% służącą do rozwiązania równań różniczkowych opisujących obiekt.
global plant_C;      plant_C = 0.75;
global plant_alpha;  plant_alpha = 15.9;

% Definicja punktu pracy. 
% Podobnie jak w przypadku definicji stałych obiektu, deklarowane są jako
% zmienne globalne.
global plant_T_C0;   plant_T_C0 = 16.97;
global plant_T_H0;   plant_T_H0 = 74.41;
global plant_T_D0;   plant_T_D0 = 35.31;
global plant_F_C0;   plant_F_C0 = 32;
global plant_F_H0;   plant_F_H0 = 19;
global plant_F_D0;   plant_F_D0 = 7;
global plant_tau_C0; plant_tau_C0 = 100;
global plant_tau0;   plant_tau0 = 55;
global plant_h0;     plant_h0 = 13.3;
global plant_T0;     plant_T0 = 38;
global plant_V0;     plant_V0 = plant_C*plant_h0^3;
