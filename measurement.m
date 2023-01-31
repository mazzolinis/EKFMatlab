sigma_w = 8e-4*sqrt(dt); % rumore additivo sul giroscopio 
sigma_f = 1e-1*sqrt(dt); % rumore additivo sull'accelerometro
sigma_bw = 1e-2*sqrt(dt); % random walk del bias del giroscopio
sigma_bf = 1e-1*sqrt(dt); % random walk del bias dell'accelerometro
%% prova senza rumore 
% sigma_w = 0;
% sigma_f = 0;
% sigma_bw = 0;
% sigma_bf = 0; 
%%
sigma_a = 3e-2*sqrt(dt); % rumore degli encoder
sigma_s = 1e-2*sqrt(dt); % rumore di posizione delle gambe
g = [0 0 -9.81]';

f_meas = NaN(3,k_max + 1);
alpha = NaN(3,4,k_max + 1);
w_meas = NaN(3,k_max + 1);
bw_real = NaN(3,k_max + 2);
bw_real(:,1) = [0 0 0]';
bf_real = bw_real;
s_meas = s_real;



for k = 0:k_max
    % valore reale dell'accelerometro
    f_meas(:,k+1) = C(:,:,k+1)*(a_real(:,k+1) - g);

    % aggiungo il rumore
    % il termine k+1 è al passo attuale, k+2 al successivo
    bf_real(:,k+2) = bf_real(:,k+1) + sigma_bf*randn(3,1);
    f_meas(:,k+1) = f_meas(:,k+1) + bf_real(:,k+1) + sigma_f/sqrt(dt)*randn(3,1);
    bw_real(:,k+2) = bw_real(:,k+1) + sigma_bw*randn(3,1); 
    w_meas(:,k+1) = w_real(:,k+1) + bw_real(:,k+1) + sigma_w/sqrt(dt)*randn(3,1); 

    for leg = 1:4
        alpha(:,leg,k+1) = invkin(leg,s_real(:,leg,k+1),sigma_a,dt);
        s_meas(:,leg,k+1) = lkin(leg,alpha(:,leg,k+1));
    end
    
end
% uno = deg2rad([0 0 90]);
% s1 = lkin(1,uno);
% alpha1 = invkin(1,s1,sigma_a,dt); 

%% ------------------------------- functions ------------------------------

function theta = invkin(leg,s,sigma_a,dt)


% parameters (need to match with lkin)
x0 = 0.25; % half length of robot back
y0 = 0.15; % half width
a2 = 0.4;
a3 = 0.5;

if leg == 1
    s = s - [x0 y0 0]';
    delta = -pi/2;
    
elseif leg == 2
    s = s - [x0 -y0 0]';
    delta = pi/2;

elseif leg == 3
    s = s - [-x0 -y0 0]';
    delta = pi/2;

elseif leg == 4
    s = s - [-x0 y0 0]';
    delta = -pi/2;
end

R = [cos(delta) -sin(delta) 0;
        sin(delta) cos(delta) 0;
        0 0 1];
s = R*s;
px = s(1);
py = s(2);
pz = s(3);

% applico la soluzione inversa del manipolatore antropomorfo dal libro di Siciliano
t1 = atan2(py,px);
c3 = (px^2 + py^2 + pz^2 - a2^2 - a3^2)/(2*a2*a3);
s3 =  sqrt(1 - c3^2); %
t3 = atan2(s3,c3);
s2 = (pz*(a2+a3*c3)-a3*s3*(px^2+py^2))/(a2^2+a3^2+2*a2*a3*c3);
c2 = (sqrt(px^2+py^2)*(a2+a3*c3) + pz*a3*s3)/(a2^2+a3^2+2*a2*a3*c3);
t2 = atan2(s2,c2);

theta = [t1 t2 t3]' + sigma_a/sqrt(dt)*randn(3,1);
end

