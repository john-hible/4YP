%% Barrier Tester on Turbine System! (using Sondergaard et al- 2012)
%Model based off of OpenFAST 5MW turbine

%% Simplified Version to test it out
clear all
close all
echo on
%Move on to 3 state system! (set pitch)
%% Drive Train Sub system
%% States - 
% omega_r - Rotor angular velocity [rad/s]
% omega_g - Generator angular velocity [rad/s]
% feta_delta - Drive train torsional angle [rad]

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
v_wind = 15;
% Also ignoring tower top motion and lead lag bending for simplicity
% Values from paper are designed to resemble NREL 5-MW turbine from FAST

lambda_r = 2.025; % Max rotor velocity [rad/s]
lambda_delta = 441.42e-3; % Ultimate load limit of drive train torsion [rad]

%% MAKE SURE SOSTOOLs & SeDuMi (AND YALMIP?) ARE IN PATH

options.solver='SeDuMi';

syms x1 x2 x3 ;

x = [x1;x2;x3];
%Create SOS program
prog = sosprogram(x);

% Define set of monomials for input to sos poly variables - INCLUDE 0th
% Order!!!
VEC_4 = monomials(x,0:2); % B up to order 4!
VEC_2 = monomials(x,0:1); %Sigmas up to order 2

%Define initial and unsafe spaces
% Initial state
%g_x0 = x'*Q*x;

%[prog, g_x0] = sospolyvar(prog,monomials(x,2),'wscoeff');

g_x0 = x1^2/(1^2) + x2^2/(97^2) + x3^2/(0.2^2) ; % = 1
%Want r_x1 = 1, r_x2 = 97?, r_x3 = 0.2?

%Unsafe Space
g_xu_1 = -x1^2 + lambda_r^2; % <0
% Add g_xu_4 to this one ^^^

% Create family of contraints that are in union
g_xu_2 = -x3^2 + lambda_delta^2; %<0
%g_xu_3 = -(x1 - 0.5); %<0
%g_xu_4 = -(x2 - 0.5*97); %<0

xlim_1 = x1^2-5^2; %<0
xlim_2 = x2^2-(5*97)^2;%<0
xlim_3 = x3^2 - 0.4^2;%<0

%Approximated C_q coefficients
C_q_Fit = [-0.0240889231402065,0.00108816365174505,0.0234530529270319,-9.14549613647725e-06,-0.000613719856304388,-0.00158790611519949];

% Simple one
C_q_Fit = [-0.01 0 0 0 0 0];
Pitch = 40;
Beta = 90-Pitch;

%Construct the vector field (x' = f(x))
f = [(1/J_r)*(0.5*rho*A*R*v_wind^2*(C_q_Fit(1) + C_q_Fit(2).*(90-Beta) + C_q_Fit(3).*(R*x1/v_wind) + C_q_Fit(4).*(90-Beta).^2 +C_q_Fit(5).*(90-Beta).*(R*x1/v_wind) + C_q_Fit(6).*(R*x1/v_wind).^2)- B_r*x1 - K_a*x3 - B_a*(x1 - N*x2)) ;
    (1/J_g)*(K_a*N*x3 + B_a*N*(x1 - N*x2) - B_g*x2); 
    x1 - N*x2];


%Create Barrier Function B, and sigma_1 and sigma_2
[prog, B] = sospolyvar(prog,VEC_4);

[prog,sig_1] = sossosvar(prog,VEC_2);
[prog,sig_2] = sossosvar(prog,VEC_2);
[prog,sig_3] = sossosvar(prog,VEC_2);
[prog,sig_4] = sossosvar(prog,VEC_2);
[prog,sig_5] = sossosvar(prog,VEC_2);
[prog,sig_6] = sossosvar(prog,VEC_2);
[prog,sig_7] = sossosvar(prog,VEC_2);
[prog,sig_8] = sossosvar(prog,VEC_2);
[prog,sig_9] = sossosvar(prog,VEC_2);
[prog,sig_10] = sossosvar(prog,VEC_2);


%Constrain so B is positive/negative in respective regions and gradB is negative 
% Inequalities are >= 0 
prog = sosineq(prog,-B -0.1 + (g_x0-1) );
%+ sig_2*xlim_1 + sig_3*xlim_2 + sig_4*xlim_3);

prog = sosineq(prog, B -0.1 + sig_1*g_xu_1 );
%+ sig_5*xlim_1 + sig_6*xlim_2 + sig_7*xlim_3);
%+ sig_5*g_xu_4);
%prog = sosineq(prog, B -0.1 + sig_2*g_xu_2) ;
%+ sig_3*g_xu_3 + sig_4*g_xu_4);

prog = sosineq(prog, - (diff(B,x1)*f(1) + diff(B,x2)*f(2) + diff(B,x3)*f(3)) );
%+ sig_8*xlim_1 + sig_9*xlim_2 + sig_10*xlim_3 );

% 
% % Constrain variables to within accurate region of coeff
% prog = sosineq(prog,x4-85); % pitch between -5 to 5 (85 to 95)
% prog = sosineq(prog,95-x4);
% 
% prog = sosineq(prog,x1-(0.4*v_wind/R)); % TSR between 4 and 22
% prog = sosineq(prog,(22*v_wind/R)-x1);


% %Define Objective Function
%prog = sossetobj(prog,coeff_1+coeff_3+coeff_6); %Minimise Tr(Q) = q1 + q3

%Solve!! (solver already called at top)
prog = sossolve(prog,options);

%Get Solution!
SOLV = sosgetsol(prog,B)
g_x = sosgetsol(prog,g_x0)

%Plot Solution

%% Barrier Function Plot 1

xx1 = -2:0.05:3;
xx2 = -500:10:500;
xx3 = -0.5:0.01:0.5;
[x1,x2] = meshgrid(xx1,xx2);

x3 = 0;

% B = - 3.99e-6.*x1.^4 - 2.669e-5.*x1.^3.*x2 - 2.323e-5.*x1.^3.*x3 - 1.327e-9.*x1.^3.*x4 - 2.678e-11.*x1.^3.*x5 - 0.0001425.*x1.^3 - 2.446e-5.*x1.^2.*x2.^2 - 6.053e-5.*x1.^2.*x2.*x3 - 8.194e-11.*x1.^2.*x2.*x4 + 6.918e-11.*x1.^2.*x2.*x5 + 3.677e-5.*x1.^2.*x2 - 7.198e-5.*x1.^2.*x3.^2 + 1.129e-9.*x1.^2.*x3.*x4 - 1.298e-9.*x1.^2.*x3.*x5 + 0.0001955.*x1.^2.*x3 - 0.006393.*x1.^2.*x4.^2 - 0.002274.*x1.^2.*x4.*x5 + 7.339e-9.*x1.^2.*x4 - 0.0068.*x1.^2.*x5.^2 + 2.12e-8.*x1.^2.*x5 + 8.591.*x1.^2 - 1.455e-5.*x1.*x2.^3 - 4.656e-5.*x1.*x2.^2.*x3 - 2.062e-10.*x1.*x2.^2.*x4 + 6.44e-11.*x1.*x2.^2.*x5 + 3.305e-5.*x1.*x2.^2 - 6.181e-5.*x1.*x2.*x3.^2 - 3.956e-10.*x1.*x2.*x3.*x4 + 1.704e-10.*x1.*x2.*x3.*x5 + 2.87e-6.*x1.*x2.*x3 - 0.0001249.*x1.*x2.*x4.^2 - 4.374e-5.*x1.*x2.*x4.*x5 + 3.595e-10.*x1.*x2.*x4 - 0.0001711.*x1.*x2.*x5.^2 + 4.418e-13.*x1.*x2.*x5 + 0.01497.*x1.*x2 - 9.223e-5.*x1.*x3.^3 + 5.123e-9.*x1.*x3.^2.*x4 - 3.362e-10.*x1.*x3.^2.*x5 + 0.0003412.*x1.*x3.^2 - 0.006026.*x1.*x3.*x4.^2 + 0.002485.*x1.*x3.*x4.*x5 - 8.501e-9.*x1.*x3.*x4 - 0.006914.*x1.*x3.*x5.^2 - 1.146e-8.*x1.*x3.*x5 + 1.595.*x1.*x3 + 3.476e-9.*x1.*x4.^3 + 2.163e-10.*x1.*x4.^2.*x5 + 0.04681.*x1.*x4.^2 - 3.945e-9.*x1.*x4.*x5.^2 + 0.007792.*x1.*x4.*x5 + 2.048e-7.*x1.*x4 - 7.323e-10.*x1.*x5.^3 + 0.0871.*x1.*x5.^2 + 3.395e-7.*x1.*x5 + 0.001009.*x1 + 1.003e-7.*x2.^4 - 1.799e-5.*x2.^3.*x3 + 3.462e-12.*x2.^3.*x4 + 1.343e-12.*x2.^3.*x5 - 1.863e-7.*x2.^3 - 2.033e-5.*x2.^2.*x3.^2 - 2.542e-10.*x2.^2.*x3.*x4 + 7.177e-11.*x2.^2.*x3.*x5 - 1.775e-5.*x2.^2.*x3 - 2.682e-9.*x2.^2.*x4.^2 + 4.003e-7.*x2.^2.*x4.*x5 - 3.181e-11.*x2.^2.*x4 - 3.623e-7.*x2.^2.*x5.^2 + 6.969e-11.*x2.^2.*x5 + 0.0001789.*x2.^2 - 3.357e-5.*x2.*x3.^3 - 4.114e-10.*x2.*x3.^2.*x4 + 1.163e-10.*x2.*x3.^2.*x5 - 1.022e-5.*x2.*x3.^2 - 0.0001422.*x2.*x3.*x4.^2 - 4.876e-5.*x2.*x3.*x4.*x5 - 3.27e-10.*x2.*x3.*x4 - 0.0001855.*x2.*x3.*x5.^2 - 6.039e-11.*x2.*x3.*x5 - 0.01244.*x2.*x3 + 8.071e-12.*x2.*x4.^3 - 2.153e-11.*x2.*x4.^2.*x5 + 0.0001193.*x2.*x4.^2 - 2.044e-11.*x2.*x4.*x5.^2 - 0.000101.*x2.*x4.*x5 - 1.22e-10.*x2.*x4 + 2.523e-12.*x2.*x5.^3 + 0.0004753.*x2.*x5.^2 + 6.378e-10.*x2.*x5 + 3.376e-6.*x2 - 3.154e-5.*x3.^4 + 1.369e-9.*x3.^3.*x4 - 5.38e-10.*x3.^3.*x5 + 0.0002092.*x3.^3 - 0.002285.*x3.^2.*x4.^2 + 0.0008615.*x3.^2.*x4.*x5 - 5.881e-9.*x3.^2.*x4 - 0.002881.*x3.^2.*x5.^2 - 8.051e-9.*x3.^2.*x5 + 9.085.*x3.^2 + 2.783e-9.*x3.*x4.^3 + 2.114e-9.*x3.*x4.^2.*x5 + 0.009463.*x3.*x4.^2 - 7.792e-10.*x3.*x4.*x5.^2 - 0.03435.*x3.*x4.*x5 + 3.55e-7.*x3.*x4 - 2.226e-10.*x3.*x5.^3 + 0.03882.*x3.*x5.^2 + 1.419e-7.*x3.*x5 - 0.0003923.*x3 - 4.093e-6.*x4.^4 - 2.988e-5.*x4.^3.*x5 + 6.185e-9.*x4.^3 - 3.718e-5.*x4.^2.*x5.^2 + 1.871e-8.*x4.^2.*x5 + 2.84.*x4.^2 - 2.224e-5.*x4.*x5.^3 + 4.32e-8.*x4.*x5.^2 + 4.243.*x4.*x5 + 9.796e-10.*x4 + 1.642e-6.*x5.^4 + 9.759e-9.*x5.^3 + 7.686.*x5.^2 + 1.439e-9.*x5 - 0.4334;
% B = - 1.307e-5.*x1.^4 - 6.724e-5.*x1.^3.*x2 - 1.52e-5.*x1.^3.*x3 - 6.309e-10.*x1.^3.*x4 - 7.749e-10.*x1.^3.*x5 + 0.0003705.*x1.^3 - 3.115e-5.*x1.^2.*x2.^2 - 0.0001144.*x1.^2.*x2.*x3 - 2.721e-11.*x1.^2.*x2.*x4 - 1.215e-10.*x1.^2.*x2.*x5 + 0.001094.*x1.^2.*x2 - 0.0001863.*x1.^2.*x3.^2 + 3.217e-9.*x1.^2.*x3.*x4 + 2.399e-9.*x1.^2.*x3.*x5 + 0.000396.*x1.^2.*x3 - 0.0006524.*x1.^2.*x4.^2 - 0.0001734.*x1.^2.*x4.*x5 + 1.41e-9.*x1.^2.*x4 - 0.0005426.*x1.^2.*x5.^2 + 2.121e-9.*x1.^2.*x5 + 0.315.*x1.^2 - 1.627e-5.*x1.*x2.^3 - 4.697e-5.*x1.*x2.^2.*x3 + 3.406e-11.*x1.*x2.^2.*x4 + 1.274e-10.*x1.*x2.^2.*x5 + 0.0005314.*x1.*x2.^2 - 0.0001031.*x1.*x2.*x3.^2 - 3.589e-11.*x1.*x2.*x3.*x4 + 5.758e-11.*x1.*x2.*x3.*x5 + 4.766e-5.*x1.*x2.*x3 - 7.691e-5.*x1.*x2.*x4.^2 - 1.383e-5.*x1.*x2.*x4.*x5 + 8.79e-10.*x1.*x2.*x4 - 9.517e-5.*x1.*x2.*x5.^2 - 5.627e-10.*x1.*x2.*x5 - 0.02175.*x1.*x2 - 0.000273.*x1.*x3.^3 + 1.834e-9.*x1.*x3.^2.*x4 + 6.506e-9.*x1.*x3.^2.*x5 - 5.105e-5.*x1.*x3.^2 - 0.0007254.*x1.*x3.*x4.^2 + 0.0002105.*x1.*x3.*x4.*x5 - 1.824e-9.*x1.*x3.*x4 - 0.0007976.*x1.*x3.*x5.^2 - 1.267e-9.*x1.*x3.*x5 + 0.08458.*x1.*x3 - 3.092e-9.*x1.*x4.^3 - 1.153e-8.*x1.*x4.^2.*x5 + 0.0003119.*x1.*x4.^2 - 3.94e-9.*x1.*x4.*x5.^2 + 9.753e-5.*x1.*x4.*x5 + 1.241e-9.*x1.*x4 + 8.207e-10.*x1.*x5.^3 + 0.0002073.*x1.*x5.^2 + 3.995e-9.*x1.*x5 + 0.000122.*x1 + 3.223e-7.*x2.^4 - 2.055e-5.*x2.^3.*x3 + 5.395e-13.*x2.^3.*x4 - 2.724e-12.*x2.^3.*x5 - 3.218e-6.*x2.^3 - 1.129e-5.*x2.^2.*x3.^2 + 1.543e-11.*x2.^2.*x3.*x4 + 1.035e-10.*x2.^2.*x3.*x5 - 0.0001683.*x2.^2.*x3 + 1.013e-7.*x2.^2.*x4.^2 + 1.911e-7.*x2.^2.*x4.*x5 - 3.078e-12.*x2.^2.*x4 - 1.721e-7.*x2.^2.*x5.^2 + 1.845e-12.*x2.^2.*x5 + 9.336e-5.*x2.^2 - 7.506e-5.*x2.*x3.^3 - 4.006e-11.*x2.*x3.^2.*x4 + 2.207e-10.*x2.*x3.^2.*x5 - 0.0001068.*x2.*x3.^2 - 9.466e-5.*x2.*x3.*x4.^2 - 1.88e-5.*x2.*x3.*x4.*x5 - 4.425e-10.*x2.*x3.*x4 - 0.0001167.*x2.*x3.*x5.^2 + 6.465e-11.*x2.*x3.*x5 - 0.003449.*x2.*x3 - 6.39e-12.*x2.*x4.^3 - 1.825e-11.*x2.*x4.^2.*x5 + 7.816e-7.*x2.*x4.^2 - 7.535e-11.*x2.*x4.*x5.^2 - 5.278e-7.*x2.*x4.*x5 + 1.864e-12.*x2.*x4 - 2.334e-12.*x2.*x5.^3 + 2.138e-6.*x2.*x5.^2 + 5.905e-12.*x2.*x5 + 3.132e-7.*x2 - 9.097e-5.*x3.^4 - 3.502e-10.*x3.^3.*x4 + 1.924e-9.*x3.^3.*x5 + 7.398e-5.*x3.^3 - 0.0002942.*x3.^2.*x4.^2 + 5.768e-5.*x3.^2.*x4.*x5 - 2.903e-10.*x3.^2.*x4 - 0.0004016.*x3.^2.*x5.^2 - 1.502e-9.*x3.^2.*x5 + 6.17.*x3.^2 + 5.017e-10.*x3.*x4.^3 - 2.907e-9.*x3.*x4.^2.*x5 + 0.0001495.*x3.*x4.^2 - 5.0e-9.*x3.*x4.*x5.^2 - 3.546e-5.*x3.*x4.*x5 + 3.232e-10.*x3.*x4 - 3.387e-10.*x3.*x5.^3 + 0.0001497.*x3.*x5.^2 - 6.905e-11.*x3.*x5 + 0.000283.*x3 + 4.102e-6.*x4.^4 - 4.053e-5.*x4.^3.*x5 + 3.541e-10.*x4.^3 - 2.666e-5.*x4.^2.*x5.^2 + 1.673e-9.*x4.^2.*x5 + 0.0003591.*x4.^2 - 1.409e-5.*x4.*x5.^3 + 7.674e-10.*x4.*x5.^2 + 0.0005135.*x4.*x5 - 3.902e-13.*x4 + 2.037e-5.*x5.^4 - 7.934e-11.*x5.^3 + 0.001284.*x5.^2 - 1.218e-10.*x5 - 1.1;
B = 0.4948.*x1.^2 - 0.01015.*x1.*x2 + 0.237.*x1.*x3 + 0.01592.*x1 + 5.243e-5.*x2.^2 - 0.002392.*x2.*x3 + 7.614e-5.*x2 + 3.086.*x3.^2 + 0.04329.*x3 - 1.115;
 
figure
mesh(x1,x2,B)
hold on

lengt = 101;
CO(:,:,1) = ones(lengt); % red
CO(:,:,2) = zeros(lengt); % green
CO(:,:,3) = zeros(lengt); % blue
mesh(zeros(lengt)+2.025,x2,x1*10,CO)
hold on
CO(:,:,1) = zeros(lengt); % red
CO(:,:,2) = ones(lengt); % green
CO(:,:,3) = zeros(lengt); % blue
mesh(x1,x2,zeros(lengt)+0.1,CO)

xlim([-3 3])
ylim([-200 200])

xlabel('Rotor Angular Velocity (rad/s)')
ylabel('Generator Angular Velocity (rad/s)')
zlabel('Barrier Function')
title('Simplified Drive Train Subsystem Barrier Function (Torsion angle = 0)')

%% Barrier function plot 2
[x1,x3] = meshgrid(xx1,xx3);

x2 = 50;
B = 0.6976.*x1.^2 - 0.005407.*x1.*x2 + 0.5095.*x1.*x3 + 0.3309.*x1 + 4.476e-5.*x2.^2 - 0.009285.*x2.*x3 + 0.001492.*x2 + 3.188.*x3.^2 - 0.1675.*x3 - 1.318;
 

figure
mesh(x1,x3,B)
hold on

lengt = 101;
CO(:,:,1) = ones(lengt); % red
CO(:,:,2) = zeros(lengt); % green
CO(:,:,3) = zeros(lengt); % blue
mesh(zeros(lengt)+2.025,x3,x1*10,CO)
hold on
CO(:,:,1) = zeros(lengt); % red
CO(:,:,2) = ones(lengt); % green
CO(:,:,3) = zeros(lengt); % blue
mesh(x1,x3,zeros(lengt)+0.1,CO)

xlim([-3 3])
ylim([-1 1])

xlabel('Rotor Angular Velocity (rad/s)')
ylabel('Torsional Angle (rad)')
zlabel('Barrier Function')
title('Simplified Drive Train Subsystem Barrier Function (Generator Speed = 50?)')



figure
x3 = 0;

% REMEMBER THE '-1'IN g_x0
g_x0 = @(x1,x2) -1 + 0.3159*x1.^2 - 0.02004*x1.*x2 + 0.07848*x1.*x3 + 1.657e-8*x1.*Beta + 1.497e-8*x1.*x5 + 0.00122*x2.^2 - 0.003586*x2.*x3 - 1.125e-9*x2.*Beta - 1.178e-9*x2.*x5 + 6.17*x3.^2 - 1.673e-9*x3.*Beta - 1.374e-9*x3.*x5 + 0.0008966*Beta.^2 + 0.00079*Beta.*x5 + 0.002201*x5.^2;

% New g_x0
%g_x0 = @(x1,x2) -1 - 5.028.*x1.^2 - 0.1023.*x1.*x2 - 3.544.*x1.*x3 - 3.183e-9.*x1.*x4 + 2.301e-8.*x1.*x5 + 0.8543.*x2.^2 - 0.06723.*x2.*x3 - 1.241e-11.*x2.*x4 + 2.157e-10.*x2.*x5 - 0.03007.*x3.^2 + 4.977e-9.*x3.*x4 + 1.707e-8.*x3.*x5 + 1.118.*x4.^2 + 0.6253.*x4.*x5 + 2.086.*x5.^2;
g_x0 = @(x1,x2,x3) -1 + x1.^2 + 0.0001063*x2.^2 + 25.0*x3.^2;
 
fimplicit3(g_x0,[-2.5 2.5 -120 120 -0.5 0.5])
%[-4 4, -4000 4000, -90 90]);
% 


%Plot Solution
xx1 = -2:0.05:3;
xx2 = -500:10:500;
xx3 = -0.5:0.01:0.5;
[x1,x2] = meshgrid(xx1,xx2);

hold on
lengt = 101;
CO(:,:,1) = ones(lengt); % red
CO(:,:,2) = zeros(lengt); % green
CO(:,:,3) = zeros(lengt); % blue
mesh(zeros(lengt)+2.025,x2,x1*8,CO)
hold on
CO(:,:,1) = ones(lengt); % red
CO(:,:,2) = zeros(lengt); % green
CO(:,:,3) = zeros(lengt); % blue
mesh(x1,x2,zeros(lengt)+0.44142,CO)
ylim([-120 120])
zlim([-0.5 0.5])


xlabel('Rotor Angular Velocity (rad/s)')
ylabel('Generator Angular Velocity (rad/s)')
zlabel('Drive Train Torsion angle (rad)')
title('Simplified Drive Train Subsystem Operating Region')



% Q = [0.23, -4e-8, -3.1e-8, 7.38e-9 -3.80e-8;
% -4e-8 5.5e-7 -1.16e-10 -1.30e-11 1.07e-12;
% -3.1e-8 -1.16e-10 861.7 -3.41e-12 5.97e-13;
% 7.38e-9 -1.30e-11 -3.41e-12 0.56e-4 -3.23e-11;
% -3.80e-8 1.07e-12 5.97e-13 -3.23e-11 0.13e-3];
% 
% figure
% x3 = 0; % Pitch = 10 so, feta = (90 - pitch) = 80
% x5 = 0;
% % REMEMBER THE '-1'IN g_x0
% g_x0 = @(x1,x2,x4) -1 + [x1, x2, x3, x4,x5]*Q*[x1, x2, x3,x4,x5]';
% fimplicit3(g_x0,[-2.5 2.5 -2000 2000 -120 120])



