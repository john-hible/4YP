function dydt = sys_norm1(t,y)

%% Drive Train Sub system
% States - 
% omega_r - Rotor angular velocity [rad/s]
% omega_g - Generator angular velocity [rad/s]
% feta_delta - Drive train torsional angle [rad]

%feta_fa - tower top foreaft bending angle (in direction of wind)
%omega_fa - Tower top foreaft bending angular velocity (in direction of wind)

% feta_beta - Blade-pitch [deg] (actually 90 minus pitch)

%% Constants
J_r = 19.38e6; % Rotor inertia [kg m2]
B_r = 150e3; % Rotor friction [Nm/(rad/s)]
K_a = 15e6; % Drive train total spring constant [Nm/rad]
B_a = 12.15e6; % Drive train total torsional friction [Nm/(rad/s)]
N = 1/97; % Gear ratio
J_g = 534.12; % Generator inertia [kg m2]
B_g = 0; % Generator friction [Nm/(rad/s)] (Frictionless)

%Tower Top
K_fa = 10.13e9; % spring constant of the hinge 
B_fa = 221e6; % damping in the hinge
L = 87.6; %Tower Length
M_n = 557.97e3; %Mass of tower top

P_rated = 5.29661e6;
Tau_rated = 43094;

%PI controller
% KP(? = 0º) = 0.01882681 s, KI(? = 0º) = 0.008068634 (from definition of
% 5MW) 
ref = 12.1*2*pi/60;
Kp = 50;
Ki = 5;
% Kp = 0.01882681;
% Ki = 0.008068634;

%Ignore Wind Disturbance for now! (otherwise look up in table)
rho = 1.22521;
A = 11.88e3;
R = 61.5;
F_g = M_n*9.81;

v_wind = 15;
v_w_eff = v_wind - L*y(5);

% Also ignoring tower top motion and lead lag bending for simplicity
% Values from paper are designed to resemble NREL 5-MW turbine from FAST
C_q_Fit = [-0.00845641270737012,0.000687609914773001,0.0179412949214744,-7.81146089582775e-06,-0.000708685201549979,-0.00111642539872539];
C_t_Fit = [-0.271077037770784,0.0116800495619584,0.197049966000335,-0.000105359159761962,-0.00771319422325017,-0.00836914583446118];

C_q_Fit = [-0.0240889231402065,0.00108816365174505,0.0234530529270319,-9.14549613647725e-06,-0.000613719856304388,-0.00158790611519949];
C_t_Fit = [-0.407703867241420,0.0153920094483093,0.244017959995536,-0.000121004122118517,-0.00687654914812738,-0.0124188251797834];
%


TSR = R*y(1)/v_w_eff;
Pitch = 90 - y(6);

C_q = C_q_Fit(1) + C_q_Fit(2).*Pitch + C_q_Fit(3).*TSR + C_q_Fit(4).*Pitch.^2 +C_q_Fit(5).*Pitch.*TSR + C_q_Fit(6).*TSR.^2;
C_t = C_t_Fit(1) + C_t_Fit(2).*Pitch + C_t_Fit(3).*TSR + C_t_Fit(4).*Pitch.^2 +C_t_Fit(5).*Pitch.*TSR + C_t_Fit(6).*TSR.^2;

F_aero =0*0.5*rho*A*v_w_eff^2*C_t;
Tau_aero = 0.5*rho*A*R*v_w_eff^2*C_q;

omega_r_dot = (1/J_r)*(-4.1827e6 + Tau_aero - B_r*y(1) - 0*K_a*y(3) - 0*B_a*(y(1) - N*y(2)));

feta_beta_dot = -Kp*omega_r_dot-Ki*(y(1)-ref);

% if feta_beta_dot > 8 
%     feta_beta_dot = 8;
% end
% 
% if feta_beta_dot < -8 
%     feta_beta_dot = -8;
% end


dydt = [omega_r_dot;
    0*(1/J_g)*(K_a*N*y(3) + B_a*N*(y(1) - N*y(2))-P_rated/y(2)); 
    0*y(1) - 0*N*y(2);
    y(5);
    (1/(M_n*L))*(-B_fa*y(5)/L + ((F_g*L-K_fa)/L)*y(4) + F_aero);
    feta_beta_dot];




