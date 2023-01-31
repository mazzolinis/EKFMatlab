clear all
close all
clc

t_max = 20;
dt = 0.01;
lx = 0.1;
% lx = 0;
f_x = 0.25;
f_y = 0;
ly = 0;
theta = deg2rad(20);
% theta = 0;
f_wx = 0.05;
omega = 2*pi*f_x;

real_movement; % sistema reale

figure(1)
plot(q_real')
title('Quaternioni')

figure(2)
plot(r_real')
title('Posizione reale del robot')
legend('x', 'y', 'z')

%% --------------------------- ERRORI DI MISURA ---------------------------

close all
measurement;

figure(1)
plot(squeeze(s_real(:,1,:))')
hold on
plot(squeeze(s_meas(:,1,:))')

figure(2)
plot(w_meas')
legend('wx','wy','wz')

figure(3)
plot(f_meas')
legend('fx','fy','fz')

figure(4)
plot(squeeze(alpha(:,1,:))');
legend('alpha 1', 'alpha 2', 'alpha 3')

%% --------------------------- STIMATORE ----------------------------------
EKF;
% trackEKF;
% navEKF;

close all
figure(1)
plot(q_error')
title('Errore di assetto')
figure(2)
plot(r_real' - r_hat')
title('Errore di posizione')
legend('x', 'y', 'z')
figure(3)
plot(q_real')
hold on
plot(q_hat')
