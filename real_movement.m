k_max = floor(t_max/dt);

r_real = zeros(3,k_max + 1);
v_real = zeros(3,k_max + 1);
a_real = zeros(3,k_max + 1);
q_real = zeros(4,k_max + 1);
w_real = zeros(3,k_max + 1);
C = zeros(3,3,k_max + 1);
p = zeros(3,4); % piedi fermi in contatto con il terreno 
s_real = zeros(3,4,k_max + 1); 

for leg = 1:4
    p(:,leg) = lkin(leg,[0 -pi/6 -pi/3]);
end

[r_real, v_real, a_real] = avanti_indietro(lx,k_max,f_x,r_real,v_real,a_real,dt);
[r_real, v_real, a_real] = destra_sinistra(ly,k_max,f_y,r_real,v_real,a_real,dt);
[q_real, w_real] = ondeggio_x(theta,k_max,f_wx,q_real,w_real,dt);
% [q_real, w_real] = rotazione_costante(k_max,f_wx,q_real,w_real,dt);

for k = 0:k_max
    C(:,:,k+1) = rotation(q_real(:,k+1));
    
    for leg = 1:4
        s_real(:,leg,k+1) = C(:,:,k+1)*(p(:,leg) - r_real(:,k+1));
    end
end


%% ------------------------ functions -------------------------------------

function [r, v, a] = avanti_indietro(l0,k_max,f,r,v,a,dt)

omega = 2*pi*f;
for k = 0:k_max
    r(1,k+1) = l0*sin(omega*k*dt);
    v(1,k+1) = omega*l0*cos(omega*k*dt);
    a(1,k+1) = -omega^2*l0*sin(omega*k*dt);
end
end

function [r, v, a] = destra_sinistra(l0,k_max,f,r,v,a,dt)

omega = 2*pi*f;
for k = 0:k_max
    r(2,k+1) = l0*sin(omega*k*dt);
    v(2,k+1) = omega*l0*cos(omega*k*dt);
    a(2,k+1) = -omega^2*l0*sin(omega*k*dt);
end
end

function [q, w] = ondeggio_x(theta,k_max,f,q,w,dt)

omega = 2*pi*f;
for k = 0:k_max

    if k < k_max/2
%     alpha = theta*sin(omega*k*dt); % angolo con moto armonico
    q(1,k+1) = sin(alpha/2);
    q(4,k+1) = cos(alpha/2);
    w(1,k+1) = omega*theta*cos(omega*k*dt);  
end
end

function [q, w] = rotazione_costante(k_max,f,q,w,dt)

omega = 2*pi*f;

for k = 0:k_max
    q(1,k+1) = sin(omega*k*dt);
    q(4,k+1) = cos(omega*k*dt);
    w(1,k+1) = omega;
end
end


