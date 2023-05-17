%Trajectory Plotter 



%% Define System
%Get sys for the model and then use [y,t,x] = initial(sys,x0) to plot!

%Flapwise System

%Start on 2 state system! 
% Flapewise Blade Bending
% States - Flapwise - Blade Tip Displacement (x1) & Blade tip Velocity(x2) (makes
% sense)
% Values from paper are designed to resemble NREL 5-MW turbine from FAST

K_flap = 55.25e3; % K_flap = Flapwise blade bending spring constant [N/m]
B_flap = 250.00e3; % B_flap = Flapwise blade bending damping constant [N/(m/s)]
M_flap = 66.00e3; % M_flap = Fictitious mass of the blade beyond the break point [kg]
F_aero = 0; %Ignore Wind Disturbance for now! (otherwise look up in table)
% Also ignoring tower top motion for simplicity

A = [0 1;-K_flap/M_flap -B_flap/M_flap];
B = [0;0];
C = [1 0;0 1];
D = [0;0];


%% Pitch Actuation System - This imitates the shutdown procedure of pitch rate = 8 degrees/s
%x1 = Feta_b (90-angle), x2 = Omega_b (-pitch rate)   
% 
% A = [0 1;-0.0894 -0.6];
% B = [0;0];
% C = [1 0;0 1];
% D = [0;0];
% 

%% plug into initial
sys = ss(A,B,C,D);

% Run 'initial'
x0 = [-22,96]; % Initial Condition

[y,t,x] = initial(sys,x0);

xx1 = -120:5:120;
xx2 = -120:5:120;
[x1,x2] = meshgrid(xx1,xx2);

% % Original
% x_u_1 = nsidedpoly(4, 'Center', [11.57+30 0], 'Sidelength', 60);
% x_u_2 = nsidedpoly(4, 'Center', [-11.57-30 0], 'Sidelength', 60);
% g_x0 = @(x1,x2) 0.009531*x1.^2 + 0.001985*x1.*x2 + 0.00174*x2.^2 - 1;
% B = 0.009531*x1.^2 + 0.001985*x1.*x2 + 0.00174*x2.^2 - 1.1;
% 

% Union of inequalities example
x_u_1 = nsidedpoly(4, 'Center', [11.57+55 55], 'Sidelength', 110);
x_u_2 = nsidedpoly(4, 'Center', [200 200], 'Sidelength', 1);
g_x0 = @(x1,x2) -1 + 0.008964*x1.^2 + 0.004039*x1.*x2 + 0.000561*x2.^2;


B = 0.008964*x1.^2 + 0.004039*x1.*x2 + 1.439e-6*x1 + 0.000561*x2.^2 + 3.6e-7*x2 - 1.1;
 

%% Plot Solution


% Plot Trajectories
len = length(y);
y1 = y(1:len,1);
y2 = y(1:len,2);
plot(y1,y2)
hold on
plot(y1(1),y2(1),'mo')



[U,V2] = gradient(-B); 
%contour(x1,x2,B)
hold on
quiver(x1,x2,U,V2,'b')
hold on
plot(x_u_1, 'FaceColor', 'r')
hold on
plot(x_u_2, 'FaceColor', 'r')
hold on
fimplicit(g_x0,'Color','g')
%'MeshDensity',300,
xlim([-30 30])
ylim([-100 100])
xlabel('Flapwise Tip Displacement (m)')
ylabel('Flapwise Tip Velocity (m/s)')
title('Flapwise Blade Bending State Space, Green = operating region, Red = unsafe region')



