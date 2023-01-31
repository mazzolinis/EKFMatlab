% allocazione di memoria
a_hat = NaN(3,k_max + 1);
bw_hat = NaN(3,k_max + 1);
bf_hat = NaN(3,k_max + 1);
C_hat = NaN(3,3,k_max + 1);
q_hat = NaN(4,k_max + 1);
r_hat = NaN(3,k_max + 1);
s_hat = NaN(12,k_max + 1);
v_hat = NaN(3,k_max + 1);
w_hat = NaN(3,k_max + 1);
f_hat = NaN(3,k_max + 1);
p_hat = NaN(3,4,k_max + 1);
P = NaN(22,22,k_max + 1); % matrice di covarianza dell'errore
q_error = NaN(3,k_max);

% inizializzazione
bw_hat(:,1) = [0 0 0]';
bf_hat(:,1) = [0 0 0]';
q_hat(:,1) = [0 0 0 1]';
temp = q_prod(q_real(:,1),[-q_hat(1,1), -q_hat(2,1), -q_hat(3,1), q_hat(4,1)]');
q_error(:,1) = temp(1:3);
r_hat(:,1) = [0 0 0]';
v_hat(:,1) = [0 0 0]';
w_hat(:,1) = w_meas(:,1) - bw_hat(:,1);
f_hat(:,1) = f_meas(:,1) - bf_hat(:,1);
P(:,:,1) = 0.2*eye(22);

for leg = 1:4
    p_hat(:,leg,1) = lkin(leg,alpha(:,leg,1));
end

acc = insAccelerometer;
gyro = insGyroscope;
filter = insEKF(acc, gyro, insMotionPose);
filter.StateCovariance = P(:,:,1);

for k = 1:k_max
% 
%     if k == 763
%         aspetta = 1;
%     end
    predict(filter,dt);
    fuse(filter,acc,f_meas(:,k+1),sigma_f^2);
    w_avg = (w_meas(:,k) + w_meas(:,k+1))/2;
    [x_prop, P(:,:,k+1)] = fuse(filter,gyro,w_avg,sigma_w^2);
%     correct(filter,f_hat(:,k+1),sigma_f^2);
    q_hat(:,k+1) = x_prop(1:4);
    q_hat(:,k+1) = [q_hat(2:4,k+1); q_hat(1,k+1)];
    w_hat(:,k+1) = x_prop(5:7);
    r_hat(:,k+1) = x_prop(8:10);
    v_hat(:,k+1) = x_prop(11:13);
    a_hat(:,k+1) = x_prop(14:16);

    temp = q_prod(q_real(:,k+1),[-q_hat(1,k+1), -q_hat(2,k+1), -q_hat(3,k+1), q_hat(4,k+1)]');
    q_error(:,k+1) = temp(1:3);
end