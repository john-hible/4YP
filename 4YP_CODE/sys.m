function dydt = sys(t,y)

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

%Ignore Wind Disturbance for now! (otherwise look up in table)
rho = 1.22521;
A = 11.88e3;
R = 61.5;
v_wind = 10;
Cq_pk = 0.08;
Cq_intcpt = 60;
% Also ignoring tower top motion and lead lag bending for simplicity
% Values from paper are designed to resemble NREL 5-MW turbine from FAST


if y(1)>0 && y(1)<10*v_wind/R
    
    TSR = R*y(1)/v_wind;
    Pitch = 90 - y(4);
    C_q = (-((TSR/5)-1).^2 + 1).*(1-(Pitch/Cq_intcpt))*Cq_pk;
    else
    C_q = 0;
end


dydt = [(1/J_r)*(0.5*rho*A*R*v_wind^2*C_q - B_r*y(1) - K_a*y(3) - B_a*(y(1) - N*y(2)));
    (1/J_g)*(K_a*N*y(3) + B_a*N*(y(1) - N*y(2))); 
    y(1) - N*y(2); 
    y(5);
    -0.6*y(5) - 0.0894*y(4)];




