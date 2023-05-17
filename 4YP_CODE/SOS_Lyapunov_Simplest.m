%% Lyapunov Tester on Turbine System! (using Sondergaard et al- 2012)
%Model based off of OpenFAST 5MW turbine

%% SIMPLEST! Version to test it out
%Ignore drive train system? Gen and torsion angle are fine anyway?
%Look at rotor speed and pitch alone?
% 3 state system maybe? - autonomous pitch?
%

clear all
close all
%echo on

%Move on to 3 state system!
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
VEC_2 = monomials(x,0:2); % So that B up to order 4! (
VEC_1 = monomials(x,0:1); % this creates Sigmas up to order 2 (from x'*Q*x)


%g_xu_1 = (x1-lambda_r-0.1)^2/(0.1^2) + (x2-70)^2/(100^2) + (x3)^2/(10^2) - 1;


%Approximated C_q coefficients
%C_q_Fit = [-0.0240889231402065,0.00108816365174505,0.0234530529270319,-9.14549613647725e-06,-0.000613719856304388,-0.00158790611519949];

% No aerodynamics
% C_q_Fit = [0 0 0 0 0 0];
% 
%Linear Aerodynamics
C_q_Fit = [0, 0, -0.0111240840652831, 0 0 0];
C_q_Fit = [0.0887009882823798, -0.000994478306871089, -0.0111240840652831, 0 0 0];

%Construct the vector field (x' = f(x))
f = [(1/J_r)*(0.5*rho*A*R*v_wind^2*(C_q_Fit(2).*x2 + C_q_Fit(3).*(R*x1/v_wind) + C_q_Fit(4).*(90-x2).^2 +C_q_Fit(5).*(90-x2).*(R*x1/v_wind) + C_q_Fit(6).*(R*x1/v_wind).^2)- B_r*x1 ) ;
     x3;
    -0.6*x3 - 0.0894*x2];


% =============================================
% The Lyapunov function V(x): 
[prog,V] = sospolyvar(prog,VEC_2,'wscoeff');

% =============================================
% Next, define SOSP constraints
[prog,sig_4] = sossosvar(prog,VEC_1);
[prog,sig_5] = sossosvar(prog,VEC_1);
[prog,sig_6] = sossosvar(prog,VEC_1);

x1_lim = (x1 - -0.5)*(x1 - 3); % <0 (x1 between 0 and 3)
x2_lim = (x2 - -10)*(x2 - 105); % <0 (x2 between -5 and 100) (pitch from -15 to 100)
x3_lim = (x3 - -10)*(x3 - 10); % <0  (x3 between -10 and 10?)

% Constraint 1 : V(x) - (x1^2 + x2^2 + x3^2) >= 0
prog = sosineq(prog,V-(x1^2+x2^2+x3^2));

% Constraint 2: -dV/dx*(x3^2+1)*f >= 0
expr = -(diff(V,x1)*f(1)+diff(V,x2)*f(2)+diff(V,x3)*f(3)) + sig_4*x1_lim + sig_5*x2_lim + sig_6*x3_lim;

prog = sosineq(prog,expr);

%prog = sossetobj(prog,-(coeff_1+coeff_3*100));

%Solve!! (solver already called at top)
prog = sossolve(prog,options);

%Get Solution!
SOLV = sosgetsol(prog,V)
%Plot Solution

%% Barrier Function Plot 1
V = matlabFunction(SOLV-20);

xx1 = -2:0.05:3;
xx2 = -150:3:150;
xx3 = -10:0.2:10;
[x1,x2] = meshgrid(xx1,xx2);

figure
mesh(x1,x2,V(x1,90-x2,0))
hold on


lengt = 101;
CO(:,:,1) = ones(lengt); % red
CO(:,:,2) = zeros(lengt); % green
CO(:,:,3) = zeros(lengt); % blue
mesh(zeros(lengt)+2.025,90-x2,x1*10,CO)
hold on
CO(:,:,1) = zeros(lengt); % red
CO(:,:,2) = ones(lengt); % green
CO(:,:,3) = zeros(lengt); % blue
mesh(x1,90-x2,zeros(lengt)+0.1,CO)

xlim([-3 3])
ylim([-20 100])
zlim([-2 10])

xlabel('Rotor Angular Velocity (rad/s)')
ylabel('Initial Pitch (Degrees)')
zlabel('Barrier Function')
title('Simplest Drive Train Subsystem Barrier Function (Initial Pitch Rate = 0)')


