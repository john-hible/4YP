%Operation Region Estimator!

% Firstly for the simplest model (3 state)


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
lambda_fa = 9.54e-3*1.5; % rad/s

rho = 1.22521;
A = 1.245e4;
R = 61.5;
L = 87.6; %Tower Length
v_wind = 25; %Doesn;t do anything! (really)
%Kp = 200;

%% Solve the ODEs
dt = 0.01;
tspan = 0:dt:30;
omega_r_0 = 1.45; %region 3 = 11/12 rpm = 1.15-1.26
ref = 12.1*2*pi/60;

leny = 3001;


pitch_list = 5:0.3:35;
rotor_spd_list = 0.5:0.016:2.1;
safe_shut = zeros(101);
unsafe_shut = zeros(101);

for j = 1:101
    
    for k = 1:101
     
        rotor_spd = rotor_spd_list(j);
        pitch = pitch_list(k);

        x0 = [rotor_spd rotor_spd*97 0 0/L 0 90-pitch -8];

        %Solve for trajectory
        [t,y] = ode45(@sys_shut_fault,tspan,x0);
    
        y_max = max(y(1:leny,1));
        %foreaft_max = max(y(1:leny,4));
        
        if y_max > lambda_r %|| foreaft_max < -lambda_fa
            unsafe_shut(j,k) = 1;
        else
            safe_shut(j,k) = 1;
        end
        
    end

end

%Plot the solution
[spd_plot,pitch_plot] = meshgrid(rotor_spd_list,pitch_list);

lengt = length(pitch_list);
CO_R(:,:,1) = ones(lengt); % red
CO_R(:,:,2) = zeros(lengt); % green
CO_R(:,:,3) = zeros(lengt); % blue

CO_G(:,:,1) = zeros(lengt); % red
CO_G(:,:,2) = ones(lengt); % green
CO_G(:,:,3) = zeros(lengt); % blue

mesh(spd_plot',pitch_plot',safe_shut,CO_G)
hold on
mesh(spd_plot',pitch_plot',unsafe_shut,CO_R)

safe_shut_shft = safe_shut*[zeros(101,1),[eye(100);zeros(1,100)]];
%[countx,county] = meshgrid(1:101,1:101);

linboi = safe_shut_shft == unsafe_shut;
spd_plot_t = spd_plot';
pitch_plot_t = pitch_plot';
spds = spd_plot_t(linboi)-0.0105/2;
pitchs = pitch_plot_t(linboi);

plot(spds,pitchs)
