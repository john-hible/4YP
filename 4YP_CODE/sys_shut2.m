function dydt = sys_shut2(t,y)

%% Drive Train Sub system
% States - 
% omega_r - Rotor angular velocity [rad/s]
% omega_g - Generator angular velocity [rad/s]
% feta_delta - Drive train torsional angle [rad]
% feta_beta - Blade-pitch [deg] (actually 90 minus pitch)
% omega_beta - Blade-pitch rate [deg/s] (actually minus pitch rate)

%% Constants
J_r = 19.38e6; % Rotor inertia [kg m2]
B_r = 150e3; % Rotor friction [Nm/(rad/s)]
K_a = 15e6; % Drive train total spring constant [Nm/rad]
B_a = 12.15e6; % Drive train total torsional friction [Nm/(rad/s)]
N = 1/97; % Gear ratio
J_g = 534.12; % Generator inertia [kg m2]
B_g = 0; % Generator friction [Nm/(rad/s)] (Frictionless)

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
v_wind = 20;
Cq_pk = 0.06;
Cq_intcpt = 60;
% Also ignoring tower top motion and lead lag bending for simplicity
% Values from paper are designed to resemble NREL 5-MW turbine from FAST


C_q_Fit = [-0.00845641270737012,0.000687609914773001,0.0179412949214744,-7.81146089582775e-06,-0.000708685201549979,-0.00111642539872539];
C_t_Fit = [-0.271077037770784,0.0116800495619584,0.197049966000335,-0.000105359159761962,-0.00771319422325017,-0.00836914583446118];

TSR = R*y(1)/v_wind;
Pitch = 90 - y(4);
C_q = C_q_Fit(1) + C_q_Fit(2).*Pitch + C_q_Fit(3).*TSR + C_q_Fit(4).*Pitch.^2 +C_q_Fit(5).*Pitch.*TSR + C_q_Fit(6).*TSR.^2;
C_t = C_t_Fit(1) + C_t_Fit(2).*Pitch + C_t_Fit(3).*TSR + C_t_Fit(4).*Pitch.^2 +C_t_Fit(5).*Pitch.*TSR + C_t_Fit(6).*TSR.^2;

% 
% if y(1)>0 && y(1)<10*v_wind/R
%     
%     TSR = R*y(1)/v_wind;
%     Pitch = 90 - y(4);
%     C_q = (-((TSR/5)-1).^2 + 1).*(1-(Pitch/Cq_intcpt))*Cq_pk;
%     else
%     C_q = 0;
% end



omega_r_dot = (1/J_r)*(0.5*rho*A*R*v_wind^2*C_q - B_r*y(1) - K_a*y(3) - B_a*(y(1) - N*y(2)));

feta_beta_dot = -Kp*omega_r_dot-Ki*(y(1)-ref);

if feta_beta_dot > 8 
    feta_beta_dot = 8;
end
if feta_beta_dot < -8 
    feta_beta_dot = -8;
end


dydt = [omega_r_dot;
    (1/J_g)*(K_a*N*y(3) + B_a*N*(y(1) - N*y(2))-P_rated/y(2)); 
    y(1) - N*y(2); 
    feta_beta_dot];
