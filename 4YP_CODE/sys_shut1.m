function dydt = sys_shut1(t,y)

%% Drive Train Sub system
% States - 
% omega_r - Rotor angular velocity [rad/s]
% omega_g - Generator angular velocity [rad/s]
% feta_delta - Drive train torsional angle [rad]

%feta_fa - tower top foreaft bending angle (in direction of wind)
%omega_fa - Tower top foreaft bending angular velocity (in direction of wind)

% feta_beta - Blade-pitch [deg] (actually 90 minus pitch)
% omega_beta - Blade-pitch rate [deg/s] (actually minus pitch rate)


%% Constants
J_r = 19.38e6; % Rotor inertia [kg m2]
B_r = 150e3; % Rotor friction [Nm/(rad/s)]
K_a = 15e6;%*55; % Drive train total spring constant [Nm/rad]
B_a = 12.15e6; %*0.15; % Drive train total torsional friction [Nm/(rad/s)]
N = 1/97; % Gear ratio
J_g = 534.12; % Generator inertia [kg m2]
B_g = 0; % Generator friction [Nm/(rad/s)] (Frictionless)

%Tower Top
K_fa = 10.13e9*1.2; % spring constant of the hinge 
B_fa = 221e6*0.35; % damping in the hinge
L = 87.6; %Tower Length
M_n = 557.97e3; %Mass of tower top

% M_flap = 66e3;
% K_flap = 55.25e3;
% B_flap = 250.00e3;
% 
% K_LL = 67.5e6;
% J_LL = 12.4e6;
% B_LL = 300e6;

%Ignore Wind Disturbance for now! (otherwise look up in table)
rho = 1.22521;
A = 1.245e4;
R = 61.5;
F_g = M_n*9.81;

v_wind = 25;
v_w_eff = v_wind - L*y(5);

% Also ignoring lead lag bending for simplicity
% Values from paper are designed to resemble NREL 5-MW turbine from FAST
C_q_Fit = [-0.00845641270737012,0.000687609914773001,0.0179412949214744,-7.81146089582775e-06,-0.000708685201549979,-0.00111642539872539];
%C_t_Fit = [-0.298870172594411,0.0186831045303255,0.212534047298766,-0.000191128780351677,-0.00713429256157738,-0.0103237054607148];

C_q_Fit = 2*[-0.0240889231402065,0.00108816365174505,0.0234530529270319,-9.14549613647725e-06,-0.000613719856304388,-0.00158790611519949];
C_t_Fit = 0.5*[-0.407703867241420,0.0153920094483093,0.244017959995536,-0.000121004122118517,-0.00687654914812738,-0.0124188251797834];
%C_t_Fit = [-0.333531105111448,0.0189722913931486,0.207893385576849,-0.000170297441388086,-0.00787753211141654,-0.00916740700259423];

%C_t_Fit = [0.00918804524480671,0.000782586502128278,-0.0517312182789156,-0.000166805956250679,0.00735234124804766,0.0441117150136366,5.79049353535173e-06,2.96517482111935e-05,-0.00280604951721465,-0.00477915112883495,-7.20507343258882e-08,-3.54345369856179e-06,-2.72309109021075e-05,0.000216502706458144,0.000199339950801924,3.01021710386354e-10,1.85338627464466e-08,9.22217860114289e-07,4.48175377371743e-07,-6.67681367985541e-06,-2.42788530536534e-06];

%C_q_Fit = [0.0887009882823798, -0.000994478306871089, -0.0111240840652831, 0 0 0];


TSR = R*y(1)/v_w_eff;
Pitch = 90 - y(6);

C_q = C_q_Fit(1) + C_q_Fit(2).*Pitch + C_q_Fit(3).*TSR + C_q_Fit(4).*Pitch.^2 +C_q_Fit(5).*Pitch.*TSR + C_q_Fit(6).*TSR.^2;
C_t = C_t_Fit(1) + C_t_Fit(2).*Pitch + C_t_Fit(3).*TSR + C_t_Fit(4).*Pitch.^2 +C_t_Fit(5).*Pitch.*TSR + C_t_Fit(6).*TSR.^2;
%C_t = polyval2(C_t_Fit,Pitch,TSR);


F_aero = 0*0.5*rho*A*v_w_eff^2*C_t;
Tau_aero = 0.5*rho*A*R*v_w_eff^2*C_q;

omega_r_dot = (1/J_r)*(Tau_aero - B_r*y(1) - 0*K_a*y(3) - 0*B_a*(y(1) - N*y(2)));


dydt = [omega_r_dot;
    (1/J_g)*(K_a*N*y(3) + B_a*N*(y(1) - N*y(2)));
    y(1) - N*y(2);
    y(5);
    (1/(M_n*L))*(-B_fa*y(5)/L + ((F_g*L-K_fa)/L)*y(4) + F_aero);
    y(7);
    -0.289*y(7) - 0.0289*y(6);];
 

% dydt = [omega_r_dot;
%     (1/J_g)*(K_a*N*y(3) + B_a*N*(y(1) - N*y(2)));
%     y(1) - N*y(2);
%     y(5);
%     (1/(M_n*L))*(-B_fa*y(5)/L + ((F_g*L-K_fa)/L)*y(4) + 0.5*F_aero + K_flap*y(8));
%     y(7);
%     -0.289*y(7) - 0.0289*y(6);
%     y(9);
%     (1/M_flap)*(0.44*F_aero - K_flap*y(8)*(M_n+M_flap)/M_n - B_flap*y(9)) + B_fa*y(5)/(M_n*L) - ((F_g*L-K_fa)/(M_n*L))*y(4);
%     ((B_r+B_a)*y(1) - B_a*N*y(2) + K_a*y(3))/J_r - (0.1965*Tau_aero + B_LL*y(10) + K_LL*y(11)*(J_r+J_LL)/J_r)/J_LL;
%     y(10)];

 
% -8;
%     0];









