%% Non-Linear Trajectory Plotter

clear all
close all
%% Define System
%Get sys for the model and then use [y,t,x] = initial(sys,x0) to plot!


%% Drive Train Sub system
% States - 
% omega_r - Rotor angular velocity [rad/s]
% omega_g - Generator angular velocity [rad/s]
% feta_delta - Drive train torsional angle [rad]

%feta_fa - tower top foreaft bending angle (in direction of wind)
%omega_fa - Tower top foreaft bending angular velocity (in direction of wind)

% feta_beta - Blade-pitch [deg] (actually 90 minus pitch)
% omega_beta - Blade-pitch rate [deg/s] (actually minus pitch rate)

lambda_r = 2.025; % Max rotor velocity [rad/s]
lambda_delta = 441.42e-3; % Ultimate load limit of drive train torsion [rad]
lambda_fa = 9.54e-3; % rad/s


R = 61.5;
v_wind = 15;
L = 87.6;
%Kp = 200;

%% Simulate!
dt = 0.01;
tspan = 0:dt:30;
omega_r_0 = 1; %region 3 = 11/12 rpm = 1.15-1.26
ref = 12.1*2*pi/60;
pitch_initial = 15;

x0 = [omega_r_0 omega_r_0*97 0 0.000 0 90-pitch_initial];

t = zeros(6001,1);
y = zeros(6001,6);
v_wind_1 = 15;
v_wind_2 = 20;

%Solve first half of trajectory
[t,y] = ode45(@sys_norm1,tspan,x0);

leny_1 = length(y);

%Solve for new wind speed
% [t(3001:6001),y(3001:6001,1:6)] = ode45(@sys_norm2,30:dt:60,y(leny_1,1:6));
% 
 leny = length(y);

% plot(y(1:leny,1),y(1:leny,2))

figure
plot(t,y(1:leny,1))
hold on
%plot(t,zeros(1,leny)+441.42e-3)
xlabel('Time (s)')
ylabel('Rotor Speed (rad/s)')
title('Simulation with change in wind speed')
ylim([0 1.8])


figure
plot(t,90-y(1:leny,6))
hold on

%Plot TSR (including wind change)
v_wind_plot = [zeros(1,3001)+v_wind_1,zeros(1,3000)+v_wind_2]';
v_wind_plot = [v_wind_1];
v_w_eff = v_wind_plot - L*y(1:leny,5);

TSR_plot =  R*y(1:leny,1)./v_w_eff;

plot(t,TSR_plot)
%plot(t,zeros(1,leny)+441.42e-3)
xlabel('Time (s)')
ylabel('Pitch (degrees)')
title('Simulation with change in wind speed')
ylim([0 15])

figure
plot(t,y(1:leny,3))
hold on
plot(t,zeros(1,leny)+441.42e-3)
xlabel('Time (s)')
ylabel('Torsion Angle')
title('Simulation with change in wind speed')
ylim([0 0.5])

%Plot tower top displacement
figure
plot(t,y(1:leny,4)*L)
hold on
plot(t,zeros(1,leny)+lambda_fa*L)
xlabel('Time (s)')
ylabel('Tower Top Foreaft Displacement (m)')
title('Tower Top Displacement')


%Plot tower top velocity
figure
plot(t,y(1:leny,5)*L)
xlabel('Time (s)')
ylabel('Tower Top Foreaft Velocity (m/s)')
title('Tower Top Velocity')

%Plot effective wind speed
figure
plot(t,v_w_eff)
xlabel('Time (s)')
ylabel('Effective  Wind Speed (m/s)')
title('Effective Wind Speed')
ylim([0 25])

%% Pitch vs TSR plot
% 
C_q_Fit = [-0.0240889231402065,0.00108816365174505,0.0234530529270319,-9.14549613647725e-06,-0.000613719856304388,-0.00158790611519949];
C_t_Fit = [-0.407703867241420,0.0153920094483093,0.244017959995536,-0.000121004122118517,-0.00687654914812738,-0.0124188251797834];
%

% Find values of TSR and Coefficient values for simulation
TSR_traj = TSR_plot;
Pitch_traj = 90 - y(1:leny,6);
Cq = C_q_Fit(1) + C_q_Fit(2).*Pitch_traj + C_q_Fit(3).*TSR_traj + C_q_Fit(4).*Pitch_traj.^2 +C_q_Fit(5).*Pitch_traj.*TSR_traj + C_q_Fit(6).*TSR_traj.^2;
Ct = C_t_Fit(1) + C_t_Fit(2).*Pitch_traj + C_t_Fit(3).*TSR_traj + C_t_Fit(4).*Pitch_traj.^2 +C_t_Fit(5).*Pitch_traj.*TSR_traj + C_t_Fit(6).*TSR_traj.^2;

% Plot mesh of Cq
TSR = (0:20);
Pitch = 0:5:100;

[TSR,Pitch] = meshgrid(TSR,Pitch);
Cq_mesh =  C_q_Fit(1) + C_q_Fit(2).*Pitch + C_q_Fit(3).*TSR + C_q_Fit(4).*Pitch.^2 +C_q_Fit(5).*Pitch.*TSR + C_q_Fit(6).*TSR.^2;
Ct_mesh =  C_t_Fit(1) + C_t_Fit(2).*Pitch + C_t_Fit(3).*TSR + C_t_Fit(4).*Pitch.^2 +C_t_Fit(5).*Pitch.*TSR + C_t_Fit(6).*TSR.^2;

figure
mesh(TSR,Pitch,Cq_mesh)
hold on
%Plot trajectory
plot3(TSR_traj,Pitch_traj,Cq+0.001,'r')
hold on
%Plot start and end with circle and cross
plot3(TSR_traj(1),Pitch_traj(1),Cq(1)+0.001,'mo')
hold on
plot3(TSR_traj(leny),Pitch_traj(leny),Cq(leny)+0.001,'rx')
hold on
%Plot boundary for TSR (based on omega r)

lengt = length(Pitch);
CO(:,:,1) = ones(lengt); % red
CO(:,:,2) = zeros(lengt); % green
CO(:,:,3) = zeros(lengt); % blue

% Plot unsafe boundary
 mesh(zeros(lengt)+2.025*R/v_wind_2,Pitch,TSR/50-0.1,CO)
% hold on

xlabel('TSR')
ylabel('Pitch (degrees)')
zlabel('Torque Coefficient (approx)')
title('Trajectory of system across mesh of Torque Coefficient')
xlim([0 12])
ylim([0 50])
zlim([-0.02 0.06])


%% Plot Ct Mesh plot

figure
mesh(TSR,Pitch,Ct_mesh)
hold on
%Plot trajectory
plot3(TSR_traj,Pitch_traj,Ct+0.001,'r')
hold on
%Plot start and end with circle and cross
plot3(TSR_traj(1),Pitch_traj(1),Ct(1)+0.001,'mo')
hold on
plot3(TSR_traj(leny),Pitch_traj(leny),Ct(leny)+0.001,'rx')
hold on

% Plot unsafe boundary
mesh(zeros(lengt)+2.025*R/v_wind,Pitch,TSR/5-1,CO)
% hold on
% v_wind = 20;
% mesh(zeros(lengt)+2.025*R/v_wind,pitch,TSR/100-0.05,CO)

xlabel('TSR')
ylabel('Pitch (°)')
zlabel('Thrust Coefficient (approx)')
title('Trajectory of system across mesh of Thrust Coefficient')
xlim([0 12])
ylim([0 90])
zlim([-1.5 0.8])



