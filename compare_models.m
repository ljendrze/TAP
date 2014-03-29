plants_coefficients;

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

den = [
   1,
   - A(2,2) - A(1,1),
   A(1,1)*A(2,2) - A(1,2)*A(2,1)
]';

num11 = [
   B(1,1)*C(1,1) + B(2,1)*C(1,2),
   A(2,1)*B(1,1)*C(1,2) - A(2,2)*B(1,1)*C(1,1) + A(1,2)*B(2,1)*C(1,1) - A(1,1)*B(2,1)*C(1,2)
]';

num12 = [
   B(1,2)*C(1,1) + B(2,2)*C(1,2),
   A(2,1)*B(1,2)*C(1,2) - A(2,2)*B(1,2)*C(1,1) + A(1,2)*B(2,2)*C(1,1) - A(1,1)*B(2,2)*C(1,2)
]';

num13 = [
   B(1,3)*C(1,1) + B(2,3)*C(1,2),
   A(2,1)*B(1,3)*C(1,2) - A(2,2)*B(1,3)*C(1,1) + A(1,2)*B(2,3)*C(1,1) - A(1,1)*B(2,3)*C(1,2)
]';

num14 = [
   B(1,4)*C(1,1) + B(2,4)*C(1,2),
   A(2,1)*B(1,4)*C(1,2) - A(2,2)*B(1,4)*C(1,1) + A(1,2)*B(2,4)*C(1,1) - A(1,1)*B(2,4)*C(1,2)
]';

num21 = [
   B(1,1)*C(2,1) + B(2,1)*C(2,2),
   A(2,1)*B(1,1)*C(2,2) - A(2,2)*B(1,1)*C(2,1) + A(1,2)*B(2,1)*C(2,1) - A(1,1)*B(2,1)*C(2,2)
]';

num22 = [
   B(1,2)*C(2,1) + B(2,2)*C(2,2),
   A(2,1)*B(1,2)*C(2,2) - A(2,2)*B(1,2)*C(2,1) + A(1,2)*B(2,2)*C(2,1) - A(1,1)*B(2,2)*C(2,2)
]';

num23 = [
   B(1,3)*C(2,1) + B(2,3)*C(2,2),
   A(2,1)*B(1,3)*C(2,2) - A(2,2)*B(1,3)*C(2,1) + A(1,2)*B(2,3)*C(2,1) - A(1,1)*B(2,3)*C(2,2)
]';

num24 = [
   B(1,4)*C(2,1) + B(2,4)*C(2,2),
   A(2,1)*B(1,4)*C(2,2) - A(2,2)*B(1,4)*C(2,1) + A(1,2)*B(2,4)*C(2,1) - A(1,1)*B(2,4)*C(2,2)
]';

tf_system = tf( {num11, num12, num13, num14; num21, num22, num23, num24}, ...
                den, ...
                'InputDelay', [ 0, plant_tau_C0, 0, 0 ], ...
                'OutputDelay', [ 0, plant_tau0 ] );

step( ss_system, tf_system );
