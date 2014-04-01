% Skrypt służący do obliczania wartości parametrów zlineryzowanego obiektu
% zbiornika z mieszaniem. Obliczane są współczynniki macierzy A, B, C, D
% modelu w postaci równań stanu, a następnie na ich podstawie formułowana
% jest postać transmitancyjna. 
% 
% Równania stanu jak i postać transmitancyjna poddawane są eksperymentowi
% skokowemu i wykreślane są wyniki tego eksperymentu w celu weryfikacji
% zgodności obydwu modeli.
% 
% W końcowej części skryptu uzyskano za pomocą polecenia c2d kilka
% transmitancji dyskretnych dla różnych okresów próbkowania. Osiągi
% wszystkich transmitancji dyskretnych porównane są na szybszym wyjściu
% obiektu, jakim jest temperatura.
% 
% Inicjalizacja zmiennych obiektu.
plants_coefficients;

% Wypełnianie macierzy równań stanu obiektu.
A = zeros(2,2);
   A(1,1) = - plant_alpha*plant_C^(-1/6)*(1/6)*plant_V0^(-5/6);
   A(1,2) = 0;
   A(2,1) = plant_V0^(-2)*( ...
               plant_F_C0*( plant_T0 - plant_T_C0 ) + ...
               plant_F_H0*( plant_T0 - plant_T_H0 ) + ...
               plant_F_D0*( plant_T0 - plant_T_D0 ) ...
            );
   A(2,2) = - ( plant_F_H0 + plant_F_C0 + plant_F_D0 )/plant_V0;

B = zeros(2,4);
   B(1,1) = 1;
   B(1,2) = 1;
   B(1,3) = 1;
   B(1,4) = 0;
   B(2,1) = ( plant_T_H0 - plant_T0 )/plant_V0;
   B(2,2) = ( plant_T_C0 - plant_T0 )/plant_V0;
   B(2,3) = ( plant_T_D0 - plant_T0 )/plant_V0;
   B(2,4) = plant_F_D0 / plant_V0;

C = zeros(2,2);
   C(1,1) = plant_C^(-1/3)*plant_V0^(-2/3)/3;
   C(1,2) = 0;
   C(2,1) = 0;
   C(2,2) = 1;

ss_system = ss(A,B,C,0,'InputDelay',[0,plant_tau_C0,0,0],'OutputDelay',[0,plant_tau0]);

% Definicja parametrów modelu w postaci transmitancyjnej. Stały mianownik
% dla wszystkich transmitancji cząstkowych i mianowniki dla poszczególnych
% transmitancji obliczane są na podstawie analitycznie wyprowadzonych
% wzorów bazujących na posiadanym modelu liniowym w postaci równań stanu.
den = [
   1,...
   - A(2,2) - A(1,1),...
   A(1,1)*A(2,2) - A(1,2)*A(2,1)
];

num11 = [
   B(1,1)*C(1,1) + B(2,1)*C(1,2),...
   A(2,1)*B(1,1)*C(1,2) - A(2,2)*B(1,1)*C(1,1) + A(1,2)*B(2,1)*C(1,1) - A(1,1)*B(2,1)*C(1,2)
];

num12 = [
   B(1,2)*C(1,1) + B(2,2)*C(1,2),...
   A(2,1)*B(1,2)*C(1,2) - A(2,2)*B(1,2)*C(1,1) + A(1,2)*B(2,2)*C(1,1) - A(1,1)*B(2,2)*C(1,2)
];

num13 = [
   B(1,3)*C(1,1) + B(2,3)*C(1,2),...
   A(2,1)*B(1,3)*C(1,2) - A(2,2)*B(1,3)*C(1,1) + A(1,2)*B(2,3)*C(1,1) - A(1,1)*B(2,3)*C(1,2)
];

num14 = [
   B(1,4)*C(1,1) + B(2,4)*C(1,2),...
   A(2,1)*B(1,4)*C(1,2) - A(2,2)*B(1,4)*C(1,1) + A(1,2)*B(2,4)*C(1,1) - A(1,1)*B(2,4)*C(1,2)
];

num21 = [
   B(1,1)*C(2,1) + B(2,1)*C(2,2),...
   A(2,1)*B(1,1)*C(2,2) - A(2,2)*B(1,1)*C(2,1) + A(1,2)*B(2,1)*C(2,1) - A(1,1)*B(2,1)*C(2,2)
];

num22 = [
   B(1,2)*C(2,1) + B(2,2)*C(2,2),...
   A(2,1)*B(1,2)*C(2,2) - A(2,2)*B(1,2)*C(2,1) + A(1,2)*B(2,2)*C(2,1) - A(1,1)*B(2,2)*C(2,2)
];

num23 = [
   B(1,3)*C(2,1) + B(2,3)*C(2,2),...
   A(2,1)*B(1,3)*C(2,2) - A(2,2)*B(1,3)*C(2,1) + A(1,2)*B(2,3)*C(2,1) - A(1,1)*B(2,3)*C(2,2)
];

num24 = [
   B(1,4)*C(2,1) + B(2,4)*C(2,2),...
   A(2,1)*B(1,4)*C(2,2) - A(2,2)*B(1,4)*C(2,1) + A(1,2)*B(2,4)*C(2,1) - A(1,1)*B(2,4)*C(2,2)
];

tf_system = tf( {num11, num12, num13, num14; num21, num22, num23, num24}, ...
                den, ...
                'InputDelay', [ 0, plant_tau_C0, 0, 0 ], ...
                'OutputDelay', [ 0, plant_tau0 ], ...
                'InputName', {'cieply strumien', 'zimny strumien', 'strumien zaklocajacy', 'temperatura zaklocajaca'}, ...
                'OutputName', {'wysokosc wody', 'temperatura wody'} );

step( ss_system, tf_system );
legend('rownania stanu', 'transmitancje','Location','East');

% Wyznaczenie czterech różnych transmitancji dyskretnych dla różnych
% okresów próbkowania celem wybrania okresu zachowującego kompromis
% pomiędzy potrzebą częstego próbkowania obiektu a dokładnością
% odwzorowania modelu ciągłego.
discrete_tf1 = c2d(tf_system,10,'zoh');
discrete_tf2 = c2d(tf_system,1,'zoh');
discrete_tf3 = c2d(tf_system,0.25,'zoh');
discrete_tf4 = c2d(tf_system,0.1,'zoh');
[Y1,T1] = step(tf_system);
[Y2,T2] = step(discrete_tf1);
[Y3,T3] = step(discrete_tf2);
[Y4,T4] = step(discrete_tf3);
[Y5,T5] = step(discrete_tf4);

figure(2);
plot(T1,Y1(:,2,1),'b');
hold on;
stairs(T2,Y2(:,2,1),'r');
stairs(T3,Y3(:,2,1),'g');
stairs(T4,Y4(:,2,1),'m');
stairs(T5,Y5(:,2,1),'c');
axis([50 250 0 0.7]);
grid on;
legend('ciagly','T_p=10','T_p=1','T_p=0.25','T_p=0.1');
xlabel('Time (seconds)');
ylabel('Amplitude');

discrete_ss = c2d(ss_system,0.25,'zoh');