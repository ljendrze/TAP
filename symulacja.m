% Skrypt słuzacy do symulacji działania obiektu jakim jest 
% zbiornik z mieszaniem. Ponadto w skrypcie tworzony jest także model
% zlinearyzowane w punkcie pracy w postaci równań stanu, a następnie
% dokonywana jest symulacja obiektu opisanego nieliniowymi równaniami
% różniczkowymi i modelu zlinearyzowanego w celu otrzymania wykresów
% obrazujących jakość linearyzacji dla podanych trajektorii zmian zmiennych
% wejściowych.
%
% Definicja stałych. 
% Są zadeklarowane jako zmienne globalne na rzecz wywoływania przez funkcję
% służącą do rozwiązania równań różniczkowych opisujących obiekt.

plants_coefficients;

% Definicja wektora czasu dla symulacji obiektu.
% Parametr dt mówi o "częstości próbkowania" ciągłego obiektu.
dt = 1;
time = [ 0 : dt : 2400 ]';

% Macierz wyjściowe obiektu, składa się z dwóch kolumn, które zawierają
% trajektorię zmian wyjść obiektu, czyli odpowiednio wysokości cieczy
% w zbiorniku i jej temperatury na wyjściu.
y_real = ones(size(time,1),2);
y_real(:,1) = plant_h0 * y_real(:,1);
y_real(:,2) = plant_T0 * y_real(:,2);

% Tablice wykorzystywane jako bufory wejściowy dla zmiennej u2 i wyjściowy
% dla zmiennej T. Bufory są wykorzystywane do opóźniania wartości zmiennych
% zgodnie z opisem działania obiektu.
u2_buffer = plant_F_C0*ones( plant_tau_C0/dt, 1 );
T_buffer = plant_T0*ones( plant_tau0/dt, 1 );

% Inicjalizacja wektora stanu początkowego (zgodnego z punktem pracy).
% Wektorem zmiennych stanu jest objętość zbiornika i temperatura mieszaniny.
plant_x0 = [ plant_C*plant_h0^3; plant_T0];

% Inicjalizacja wektora wejścia w stanie początkowym. Wejście stanowi 
% wektor temperatur wody gorącej i wody zimnej.
plant_u0 = [plant_F_H0; plant_F_C0];

% Wartość zmiennych po skoku:
% u_step = [ 19; 32 ];
u_step = [ 19; 35 ];

% Trajektoria zmian wejść obiektu. W chwili obecnej jest to stała trajektoria
% po skoku wartości zadanej po pierwszej minucie.
u_trajectory = ones(size(time,1),2);
u_trajectory(:,1) = u_trajectory(:,1)*plant_u0(1);
u_trajectory(:,2) = u_trajectory(:,2)*plant_u0(2);

u_real_trajectory = ones(size(time,1),2);
u_real_trajectory(:,1) = u_real_trajectory(:,1)*plant_u0(1);
u_real_trajectory(:,2) = u_real_trajectory(:,2)*plant_u0(2);

for i = 1 : size(time,1)
   u_trajectory(i,1) = u_step(1);
   u_real_trajectory(i,1) = u_step(1);
   u_real_trajectory(i,2) = u_step(2);

   u_trajectory(i,2) = u2_buffer(size(u2_buffer,1));
   u2_buffer = [ u_step(2); u2_buffer(1:length(u2_buffer)-1,:) ];
end



% Inicjalizacja poziomu zakłócenia, zakłócenie przyjmowane jest stałe
% dla całego przebiegu symulacji.
plant_z0 = [ plant_F_D0; plant_T_D0 ];

% Wartość zmiennych zakłócających po skokach:
% z_step = [ 7; 35.31 ];
z_step = [ 7; 35.31 ];

% Trajektoria zmian zakłóceń obiektu. W chwili obecnej jest to wartość
% stała w całym czasie trwania symulacji.
z_trajectory = ones(size(time,1),2);
z_trajectory(:,1) = z_trajectory(:,1)*plant_z0(1);
z_trajectory(:,2) = z_trajectory(:,2)*plant_z0(2);

for i = 1 : size(time,1)
   z_trajectory(i,1) = z_step(1);
   z_trajectory(i,2) = z_step(2);
end

% Opcje dodatkowe dla solwera ode45.
options = odeset('AbsTol',1e-6,'RelTol',1e-6);

u = u_trajectory(1,:)';
z = z_trajectory(1,:)';

% ======== ZMIENNE WYŚWIETLAJĄCE POSTĘP SYMULACJI ==================== %
simulation_step = 0;                                                   %
letters_written = 0;                                                   %
finish_time = size(time,1);                                            %
% ==================================================================== %

plant_V0 = plant_x0(1);

% Zlinearyzowany obiekt ciągły w postaci równań stanu:
% 
%    .
%    x = Ax + Bu
%
%    y = Cx + D
%
% Wektor stanu składa się tak jak w obiekcie, z dwóch elementów, tj.
% objętości cieczy w zbiorniku, a także jej temperatury.
% Wektor wejść składa się z trzech elementów - strumienia wody ciepłej,
% strumienia wody zimnej, a także strumienia zakłócającego.
% Zakłócenie zostało wprowadzone jako trzecie wejście.
% Wyście modelu stanowi dwuelementowy wektor - poziom cieczy w zbiorniku
% i temperatura wypływającej cieczy.
%
% WAŻNE:
% Należy pamiętać, że po linearyzacji zmienne stanu i zmienne wejściowe są
% faktycznie przyrostami wejść i stanu w odniesieniu do wartości bazowych,
% którymi są wartości odpowiadające obiektowi w punkcie pracy.
% Wobec tego wszystkie zmiany wejść przed dojściem do obiektu powinny zostać
% pomniejszone o wartości bazowe. Przykładowo: w punkcie pracy do zbiornika
% dolewana jest ciepła woda strumieniem 19cm^3/s. Zmieniając strumień
% na 22cm^3/s wykonujemy skok na wejściu o 3cm^3/s, stąd na odpowiadające
% modelowi zlinearyzowanemu wejście powinniśmy podać właśnie wartość 3.

global A;
global B;
global C;

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

D = zeros(2,1);

% Z racji na to, że model zlinearyzowany posługuje się zmianami stanu,
% wejść i wyjść w odniesieniu do punktu pracy, wyjście powinno zostać
% zainicjalizowane wartościami z punktu pracy i podczas symulacji wartość
% wyjścia z modelu obliczana z równań stanu powinna zostać dodawana
% do wyjścia reprezentowanego przez y_linearized.
y_linearized_trajectory = y_real;

sim_x0 = plant_x0;

x_linearized = [0; 0];
u_linearized = [ u(1) - plant_u0(1); 
                 u(2) - plant_u0(2);
                 z(1) - plant_z0(1);
                 z(2) - plant_z0(2)         ];
y_linearized = zeros(2,1);
x_delayed = zeros(2,1);
T_linearized_buffer = zeros( plant_tau0/dt, 1 );


for i = 2 : finish_time
   sim_t = [time(i); time(i)+dt];

   % ======== ZMIENNE WYŚWIETLAJĄCE POSTĘP SYMULACJI ================= %
   if i > simulation_step*(finish_time/100)                            %
      while letters_written > 0                                        %
         fprintf('\b');                                                %
         letters_written = letters_written - 1;                        %
      end                                                              %
      letters_written = fprintf('Progress:    %d%%',simulation_step);  %
      simulation_step = simulation_step + 1;                           %
   end                                                                 %
   % ================================================================= %

   [t, x] = ode45(@zbiornik, sim_t, sim_x0, odeset, u, z);

   sim_x0 = x(size(x,1),:);

   u = u_trajectory(i,:)';
   z = z_trajectory(i,:)';

   y_real(i,1) = (sim_x0(1)/plant_C)^(1/3);
   y_real(i,2) = T_buffer(size(T_buffer,1));

   T_buffer = [ sim_x0(2); T_buffer(1:length(T_buffer)-1,:) ];

   [t, x_linearized] = ode45( @zbiornik_linearized_state, ... 
                              sim_t,                      ...
                              x_linearized,               ...
                              options,                    ...
                              u_linearized );

   x_linearized = x_linearized(size(x_linearized,1),:);

   x_delayed(1) = x_linearized(1);
   x_delayed(2) = T_linearized_buffer(size(T_buffer,1));

   y_linearized = C * x_delayed + D;
   
   T_linearized_buffer = ...
           [ x_linearized(2); 
             T_linearized_buffer(1:length(T_linearized_buffer)-1,:) ];

   y_linearized_trajectory(i,:) = y_linearized_trajectory(i,:) + y_linearized';

   u_linearized = [ u(1) - plant_u0(1); 
                    u(2) - plant_u0(2);
                    z(1) - plant_z0(1);
                    z(2) - plant_z0(2)         ];

end

fprintf('\n\n');

if 1
   subplot(2,3,3);
   plot(time/60,y_real(:,1),'b');
   hold on;
   plot(time/60,y_linearized_trajectory(:,1),'r');
   xlabel('czas (min)');
   ylabel('wysokosc wody (cm)');
   grid on;
   legend('rzeczywisty model','model zlinearyzowany','Location','East');
   title('y_1(t)');
   
   subplot(2,3,6);
   plot(time/60,y_real(:,2),'b');
   hold on;
   plot(time/60,y_linearized_trajectory(:,2),'r');
   grid on;
   xlabel('czas (min)');
   ylabel('temperatura wody (^oC)');
   legend('rzeczywisty model','model zlinearyzowany','Location','East');
   title('y_2(t)');
   
   subplot(2,3,1);
   stairs(time/60,u_real_trajectory(:,1),'m');
   hold on;
   xlabel('czas (min)');
   ylabel('strumien goracej wody (cm^3/s)');
   grid on;
   title('u_1(t)');
   
   subplot(2,3,4);
   stairs(time/60,u_real_trajectory(:,2),'m');
   hold on;
   xlabel('czas (min)');
   ylabel('strumien zimnej wody (cm^3/s)');
   grid on;
   title('u_2(t)');
   
   subplot(2,3,2);
   stairs(time/60,z_trajectory(:,1),'m');
   hold on;
   xlabel('czas (min)');
   ylabel('strumien zaklocajacy wody (cm^3/s)');
   grid on;
   title('u_3(t)');
   
   subplot(2,3,5);
   stairs(time/60,z_trajectory(:,2),'m');
   hold on;
   xlabel('czas (min)');
   ylabel('temperatura wody zaklocajacej (^oC)');
   grid on;
   title('u_4(t)');
end
