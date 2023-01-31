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
P = NaN(27,27,k_max + 1); % matrice di covarianza dell'errore
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
P(:,:,1) = 0.2*eye(27);

for leg = 1:4
    p_hat(:,leg,1) = lkin(leg,alpha(:,leg,1));
end

sigma_p = sigma_s; % per semplicità sto tenendo tutte le zampe a terra, 
% in fase di movimento sigma_p cambierà fra le zampe a terra e quelle
% sollevate

% definisco il filtro di Kalman
% nello stato del filtro sono presenti solo le componenti vettoriali del
% quaternione (come angolo), il quarto elemento si ricava dalla norma 
stato = [r_hat(:,1); v_hat(:,1); 2*asin(q_hat(1:3,1)); 
    p_hat(:,1,1); p_hat(:,1,1); p_hat(:,1,1); p_hat(:,1,1);
    bf_hat(:,1); bw_hat(:,1)];
filtro = trackingEKF(@predizione,@misure,stato, ...
     'StateTransitionJacobianFcn',@jacobiano_predizione, ...
    'MeasurementJacobianFcn',@jacobiano_misure);
initialize(filtro, stato,P(:,:,1));

% il ciclo parte da 1 perché il primo elemento è di inizializzazione (k+1 è
% il valore attuale)
for k = 1:k_max

    C_hat(:,:,k) = rotation(q_hat(:,k));
%     w_hat(:,k+1) = w_meas(:,k+1) - bw_hat(:,k+1);
%     w_avg = (w_hat(:,k) + w_hat(:,k+1))/2; % uso la media di due misure come dal Trawny
%     f_hat(:,k+1) = f_meas(:,k+1) - bf_hat(:,k+1);
    w_hat(:,k) = w_meas(:,k) - bw_hat(:,k);
    f_hat(:,k) = f_meas(:,k) - bf_hat(:,k);
    w_avg = w_hat(:,k);

    % propagazione
    Q_k = covarianza(w_avg,C_hat(:,:,k),dt,sigma_f,sigma_bf,sigma_p,sigma_w,sigma_bw);
    filtro.ProcessNoise = Q_k;
    [x_prop, P_prop] = predict(filtro,w_avg,f_hat(:,k),C_hat(:,:,k),dt,g);
    phi = x_prop(7:9);
    q_hat(:,k+1) = [sin(phi/2); sqrt(1-norm(sin(phi/2))^2)];
% 
%     % correzione
%     [s_hat(:,k+1),R] = cinematica(alpha(:,:,k+1),sigma_s,sigma_a,dt);
%     filtro.MeasurementNoise = R;
%     [stato, P(:,:,k+1)] = correct(filtro,s_hat(:,k+1));
%     

    stato = x_prop;
    P(:,:,k+1) = P_prop;
    % correggo i quaternioni utilizzando l'errore moltiplicativo
    dphi = stato(7:9) - x_prop(7:9);
    dq = [dphi/2; 1];
    q_hat(:,k+1) = q_prod(dq,q_hat(:,k+1));
    q_hat(:,k+1) = q_hat(:,k+1)/norm(q_hat(:,k+1));
%     phi = x_corr(7:9);
%     q_hat(:,k+1) = [sin(norm(phi)/2)*phi/norm(phi); cos(norm(phi))];
    
    % aggiorno lo stato del filtro
    stato(7:9) = 2*asin(q_hat(1:3,k+1));
    filtro.State = stato;

    % aggiorno tutte le variabili di stato
    r_hat(:,k+1) = stato(1:3);
    v_hat(:,k+1) = stato(4:6);
    for leg = 1:4
        p_hat(:,leg,k+1) = stato(10+3*(leg-1):9+3*leg);
    end
    bf_hat(:,k+1) = stato(22:24);
    bw_hat(:,k+1) = stato(25:27);
    
    % calcolo il vero errore sui quaternioni
    temp = q_prod(q_real,[-q_hat(1,k+1), -q_hat(2,k+1), -q_hat(3,k+1), q_hat(4,k+1)]');
    q_error(:,k+1) = temp(1:3);

end


%% -------------------------- funzioni ------------------------------------

function x_prop = predizione(stato,w,f,C,dt,g)
bw = stato(25:27); 
bf = stato(22:24);
p = stato(10:21); % bias e posizioni rimangono invariate in predizione
r = stato(1:3);
v = stato(4:6);
phi = stato(7:9);
q = [sin(phi/2); sqrt(1-norm(sin(phi/2))^2)];
a = C'*f + g;
v_prop = v + a*dt;
r_prop = r + v*dt + a*dt^2/2;
% omega_m = omega_mat(w); 
% q_prop = expm(omega_m*dt/2)*q;
q_prop = q_prod([w*dt/2; 0], q);
phi_prop = 2*asin(q_prop(1:3));

x_prop = [r_prop; v_prop; phi_prop; p; bf; bw];
end

function F = jacobiano_predizione(stato,w,f,C,dt,g)
% calcolo di F_k
G1 = gamma(w,1,dt);
F13 = -dt^2*C'*skew(f)/2;
F23 = -dt*C'*skew(f);
F15 = -dt^2*C'/2;
F25 = -dt*C';
G0 = expm(dt*skew(dt*w));
F_up = [eye(3), dt*eye(3), F13, zeros(3,12), F15,zeros(3);
    zeros(3), eye(3), F23, zeros(3,12), F25, zeros(3);
    zeros(3), zeros(3), G0', zeros(3,12), zeros(3), -G1'];
F = [F_up; zeros(18,9), eye(18)];
end

function G = gamma(w,n,dt)
G = 0;
for i = 1:5
    G = G + dt^(i+n)*skew(w)^i/factorial(i+n);
end
end

function h = misure(stato)

phi = stato(7:9);
q = [sin(phi/2); sqrt(1-norm(sin(phi/2))^2)];
C = rotation(q);
r = stato(1:3);
p1 = stato(10:12);
p2 = stato(13:15);
p3 = stato(16:18);
p4 = stato(19:21);
h = [C*(p1-r); C*(p2-r); C*(p3-r); C*(p4-r)];
end

function H = jacobiano_misure(stato)

phi = stato(7:9);
q = [sin(phi/2); sqrt(1-norm(sin(phi/2))^2)];
C = rotation(q);
r = stato(1:3);
p1 = stato(10:12);
p2 = stato(13:15);
p3 = stato(16:18);
p4 = stato(19:21);
H = [-C, zeros(3), skew(C*(p1-r)); 
    -C, zeros(3), skew(C*(p2-r)); 
    -C, zeros(3), skew(C*(p3-r));
    -C, zeros(3), skew(C*(p4-r))];
H = [H, blkdiag(C,C,C,C)];
H = [H, zeros(12,6)];
end

function Q = covarianza(w,C,dt,sigma_f,sigma_bf,sigma_p,sigma_w,sigma_bw)

Q11 = dt^3*sigma_f^2*eye(3)/3 + dt^5*sigma_bf^2*eye(3)/20;
Q12 = dt^2*sigma_f^2*eye(3)/2 + dt^4*sigma_bf^2*eye(3)/8;
Q15 = -dt^3*C'*sigma_bf^2*eye(3)/6;
Q21 = dt^2*sigma_f^2*eye(3)/2 + dt^4*sigma_bf^3*eye(3)/8;
Q22 = dt*sigma_f^2*eye(3) + dt^3*sigma_bf^2*eye(3)/3;
Q25 = - dt^2*C'*sigma_bf^2*eye(3)/2;
G3 = gamma(w,3,dt);
Q33 = dt*sigma_w^2*eye(3) + (G3+G3')*sigma_bf^2*eye(3);
G2 = gamma(w,2,dt);
Q36 = -G2'*sigma_bf^2*eye(3);
Q44 = dt*C'*sigma_p^2*eye(3)*C;
Q44 = blkdiag(Q44,Q44,Q44,Q44); % ho 4 zampe perciò ripeto il calcolo 4 volte
Q51 = -dt^3*sigma_bf^2*eye(3)*C/6;
Q52 = -dt^2*sigma_bf^2*C/2;
Q55 = dt*sigma_bf^2*eye(3);
Q63 = -sigma_bw^2*G2;
Q66 = dt*sigma_bw^2*eye(3);

Q = [Q11, Q12, zeros(3), zeros(3,12), Q15, zeros(3);
    Q21, Q22, zeros(3), zeros(3,12), Q25, zeros(3);
    zeros(3), zeros(3), Q33, zeros(3,12), zeros(3), Q36;
    zeros(12,9), Q44, zeros(12,6);
    Q51, Q52, zeros(3), zeros(3,12), Q55, zeros(3);
    zeros(3,6), Q63, zeros(3,15), Q66];
Q = (Q + Q')/2; % cerco di ridurre i problemi di simmetria
end

function [s,R] = cinematica(alpha, sigma_s, sigma_a, dt)

s = 0;
R = 0;
for leg = 1:4
    [s_i,J] = lkin(leg,alpha(:,leg));
    s = [s; s_i];
    R_i = sigma_s^2*eye(3) + J*sigma_a^2*J';
    R = blkdiag(R,R_i);
end
s = s(2:end) + sigma_s/sqrt(dt)*randn(12,1);
R = R(2:end,2:end);
end
