% allocazione di memoria
a_hat = NaN(3,k_max + 1);
bw_hat = NaN(3,k_max + 1);
bf_hat = NaN(3,k_max + 1);
C_hat = NaN(3,3,k_max + 1);
q_hat = NaN(4,k_max + 1);
r_hat = NaN(3,k_max + 1);
s_hat = NaN(3,4,k_max + 1);
v_hat = NaN(3,k_max + 1);
w_hat = NaN(3,k_max + 1);
f_hat = NaN(3,k_max + 1);
p_hat = NaN(3,4,k_max + 1);
P = NaN(27,27,k_max + 1); % matrice di covarianza dello stato
q_error = NaN(3,k_max);

% inizializzazione
bw_hat(:,1) = [0 0 0]';
bf_hat(:,1) = [0 0 0]';
q_hat(:,1) = [0 0 0 1]';
r_hat(:,1) = [0 0 0]';
v_hat(:,1) = [0 0 0]';
w_hat(:,1) = w_meas(:,1) - bw_hat(:,1);
f_hat(:,1) = f_meas(:,1) - bf_hat(:,1);
P(:,:,1) = 0.2*eye(27);

for leg = 1:4
    p_hat(:,leg,1) = lkin(leg,alpha(:,leg,1));
end

% il ciclo parte da 1 perché il primo elemento è di inizializzazione (k+1 è
% il valore attuale)
for k = 1:k_max
    
    % predizione
    bf_hat(:,k+1) = bf_hat(:,k);
    bw_hat(:,k+1) = bw_hat(:,k);
    p_hat(:,:,k+1) = p_hat(:,:,k);
    w_hat(:,k+1) = w_meas(:,k+1) - bw_hat(:,k+1);
    f_hat(:,k+1) = f_meas(:,k+1) - bf_hat(:,k+1);
    C_hat(:,:,k) = rotation(q_hat(:,k));
    a_hat(:,k) = C_hat(:,:,k)'*(f_meas(:,k) - bf_hat(:,k)) + g;
    v_hat(:,k+1) = v_hat(:,k) + a_hat(:,k)*dt;
    r_hat(:,k+1) = r_hat(:,k) + v_hat(:,k)*dt + a_hat(:,k)*dt^2/2;

%     w_avg = (w_hat(:,k) + w_hat(:,k+1))/2;% media di due misure dal Trawny
    w_avg = w_hat(:,k+1);
    omega_m = omega_mat(w_avg); 
%     q_hat(:,k+1) = expm(omega_m*dt/2)*q_hat(:,k);
    q_hat(:,k+1) = q_prod([w_avg*dt/2;0],q_hat(:,k));

    % calcolo di F_k
    G1 = gamma(w_avg,1,dt);
    F13 = -dt^2*C_hat(:,:,k)'*skew(f_hat(:,k+1))/2;
    F23 = -dt*C_hat(:,:,k)'*skew(f_hat(:,k+1));
    F15 = -dt^2*C_hat(:,:,k)'/2;
    F25 = -dt*C_hat(:,:,k)';
    G0 = expm(dt*skew(dt*w_avg));
    F_up = [eye(3), dt*eye(3), F13, zeros(3,12), F15,zeros(3);
        zeros(3), eye(3), F23, zeros(3,12), F25, zeros(3);
        zeros(3), zeros(3), G0', zeros(3,12), zeros(3), -G1'];
    F = [F_up; zeros(18,9), eye(18)];

    % calcolo di Q_k
    Q11 = dt^3*sigma_f^2*eye(3)/3 + dt^5*sigma_bf^2*eye(3)/20;
    Q12 = dt^2*sigma_f^2*eye(3)/2 + dt^4*sigma_bf^2*eye(3)/8;
    Q15 = -dt^3*C_hat(:,:,k)'*sigma_bf^2*eye(3)/6;
    Q21 = dt^2*sigma_f^2*eye(3)/2 + dt^4*sigma_bf^3*eye(3)/8;
    Q22 = dt*sigma_f^2*eye(3) + dt^3*sigma_bf^2*eye(3)/3;
    Q25 = - dt^2*C_hat(:,:,k)'*sigma_bf^2*eye(3)/2;
    G3 = gamma(w_avg,3,dt);
    Q33 = dt*sigma_w^2*eye(3) + (G3+G3')*sigma_bf^2*eye(3);
    G2 = gamma(w_avg,2,dt);
    Q36 = -G2'*sigma_bf^2*eye(3);
    Q44 = dt*C_hat(:,:,k)'*sigma_s^2*eye(3)*C_hat(:,:,k);
    Q44 = blkdiag(Q44,Q44,Q44,Q44); % ho 4 zampe perciò ripeto il calcolo 4 volte
    Q51 = -dt^3*sigma_bf^2*eye(3)*C_hat(:,:,k)/6;
    Q52 = -dt^2*sigma_bf^2*C_hat(:,:,k)/2;
    Q55 = dt*sigma_bf^2*eye(3);
    Q63 = -sigma_bw^2*G2;
    Q66 = dt*sigma_bw^2*eye(3);

    Q = [Q11, Q12, zeros(3), zeros(3,12), Q15, zeros(3);
        Q21, Q22, zeros(3), zeros(3,12), Q25, zeros(3);
        zeros(3), zeros(3), Q33, zeros(3,12), zeros(3), Q36;
        zeros(12,9), Q44, zeros(12,6);
        Q51, Q52, zeros(3), zeros(3,12), Q55, zeros(3);
        zeros(3,6), Q63, zeros(3,15), Q66];

    P_prop = F*P(:,:,k)*F' + Q;
    P(:,:,k+1) = P_prop;

    % correzione
    y = 0;
    H = zeros(1,9);
    R = 0;
    C_hat(:,:,k+1) = rotation(q_hat(:,k+1)); % corrisponde a C_k-, verrà corretta al prossimo passo di k
    for leg = 1:4
        [s,J] = lkin(leg,alpha(:,leg,k+1));
        y = [y; s - C_hat(:,:,k+1)*(p_hat(:,leg,k+1) - r_hat(:,k+1))];
        R_i = sigma_s^2*eye(3) + J*sigma_a^2*J';
        R = blkdiag(R,R_i);
        H = [H;
            -C_hat(:,:,k+1), zeros(3), skew(C_hat(:,:,k+1)*(p_hat(:,leg,k+1) - r_hat(:,k+1)))];
    end
    
    y = y(2:end);
    H = H(2:end,:);
    R = R(2:end,2:end);
    H = [H, blkdiag(C_hat(:,:,k+1),C_hat(:,:,k+1),C_hat(:,:,k+1),C_hat(:,:,k+1))];
    H = [H, zeros(12,6)];

    S = H*P_prop*H' + R;
    K = P_prop*H'/S;
    dx = K*y;
    P(:,:,k+1) = (eye(27)-K*H)*P_prop*(eye(27)-K*H)' + K*R*K';

    % aggiorno le variabili di stato
    r_hat(:,k+1) = r_hat(:,k+1) + dx(1:3);
    v_hat(:,k+1) = v_hat(:,k+1) + dx(4:6);

    dq = [dx(7:9)/2; 1];
    q_hat(:,k+1) = q_prod(dq,q_hat(:,k+1));
    q_hat(:,k+1) = q_hat(:,k+1)/norm(q_hat(:,k+1)); % normalizzo il quaternione

    dp = [dx(10:12), dx(13:15), dx(16:18), dx(19:21)];
    p_hat(:,:,k+1) = p_hat(:,:,k+1) + dp;

    bf_hat(:,k+1) = bf_hat(:,k+1) + dx(22:24);
    f_hat(:,k+1) = f_meas(:,k+1) - bf_hat(:,k+1);
    bw_hat(:,k+1) = bw_hat(:,k+1) + dx(25:27);
    w_hat(:,k+1) = w_meas(:,k+1) - bw_hat(:,k+1);

    q_eps = q_prod(q_real(:,k), [-q_hat(1,k), -q_hat(2,k), -q_hat(3,k), q_hat(4,k)]');
    q_error(:,k) = 2*q_eps(1:3);
end


%% -------------------------------- functions -----------------------------
function G = gamma(w,n,dt)

G = 0;

for i = 0:5
    G = G + dt^(i+n)*skew(w)^i/factorial(i+n);
end

end