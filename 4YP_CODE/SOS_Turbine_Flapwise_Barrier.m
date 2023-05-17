%% Barrier Tester on Turbine System! (using Sondergaard et al- 2012)
clear all
close all

%Start on 2 state system! 
% Flapewise Blade Bending
% States - Flapwise - Blade Tip Displacement (x1) & Blade tip Velocity(x2) (makes
% sense)
% K_flap = Flapwise blade bending spring constant [N/m]
% B_flap = Flapwise blade bending damping constant [N/(m/s)]
% M_flap = Fictitious mass of the blade beyond the break point [kg]
% F_aero = Aerodynamic Thrust Force (N)
% Values from paper are designed to resemble NREL 5-MW turbine from FAST

K_flap = 55.25e3;
B_flap = 250.00e3;
M_flap = 66.00e3;
F_aero = 0; %Ignore Wind Disturbance for now! (otherwise look up in table)
% Also ignoring tower top motion for simplicity

% Max flapwise blade tip displacement = 11.57m
lambda = 11.57;

%% MAKE SURE SOSTOOLs & SeDuMi ARE IN PATH

options.solver='SeDuMi';

syms x1 x2;
x = [x1;x2];

%Create SOS program
prog = sosprogram([x1;x2]);
[prog, g_x0] = sospolyvar(prog,monomials(x,2),'wscoeff');

% pvar x1 x2;
% x = [x1;x2];
% 
% %Create SOS program
% prog = sosprogram([x1;x2]);
% [prog, g_x0] = sosquadvar(prog,x',x,1,1,'pos');

%Define initial and unsafe spaces
% Initial state
%g_x0 = x'*Q*x;

%a = 0;
%g_x0 = [x; 1]'*[q1 q2 c1/2;q2 q3 c2/2;c1/2 c2/2 d]*[x;1];

g_xu1 = -x1+lambda; % <0
g_xu2 = (x1+lambda); % <0
%g_xu3 = -x2; % <0

%Inequalities hold when this is true!
x_lim = -x2; %<0 

%Construct the vector field (x' = f(x))
f = [x2; (1/M_flap)*(F_aero - K_flap*(x1) - B_flap*x2)];

% Define set of monomials for input to sos poly variables - INCLUDE 0th
% Order!!!
VEC_4 = monomials([x1; x2],0:2); % B up to order 4!
VEC_2 = monomials([x1; x2],0:2); %Sigmas up to order 2
VEC_1 = monomials([x1; x2],0:2); % Up to order 1

%Create Barrier Function B, and sigma_1 and sigma_2
[prog, B] = sospolyvar(prog,VEC_2);

[prog,sig_1] = sossosvar(prog,VEC_1);
[prog,sig_2] = sossosvar(prog,VEC_1);

[prog,sig_3] = sossosvar(prog,VEC_1);
[prog,sig_4] = sossosvar(prog,VEC_1);
[prog,sig_5] = sossosvar(prog,VEC_1);

%Constrain so B is positive/negative in respective regions and gradB is negative 
% Inequalities are >= 0 
% If you want multiple inequalities to hold at the same time, then put them
% in the same line (with seperate sigmas)

% To limit state space

prog = sosineq(prog,-B -0.1 + (g_x0-1) + sig_3*x_lim );
%

prog = sosineq(prog, B -0.1 + sig_1*g_xu1 + sig_4*x_lim);

%prog = sosineq(prog, B -0.1 + sig_2*g_xu2 + sig_5*x_lim);

prog = sosineq(prog, - (diff(B,x1)*f(1) + diff(B,x2)*f(2)) + sig_5*x_lim);

% prog = sosineq(prog, Q1);
% prog = sosineq(prog, Q2);
% prog = sosineq(prog, Q3);


% %Define Objective Function
prog = sossetobj(prog,coeff_1+coeff_3); %Minimise Tr(Q) = q1 + q3

%Solve!! (solver already called at top)
prog = sossolve(prog,options);

%Get Solution!
SOLV = sosgetsol(prog,B)
g_x0 = sosgetsol(prog,g_x0)
obj = sosgetsol(prog,coeff_1+coeff_3)

%Plot Solution
xx1 = -30:0.5:30;
xx2 = -30:0.5:30;
[x1,x2] = meshgrid(xx1,xx2);

B = matlabFunction(SOLV);

% B = 0.009531*x1.^2 + 0.001985*x1.*x2 + 0.00174*x2.^2 - 1.1;  
% B = 0.009531*x1.^2 + 0.001985*x1.*x2 + 1.941e-11*x1 + 0.00174*x2.^2 + 1.257e-11*x2 - 1.1;
%  
figure
mesh(x1,x2,B(x1,x2))
hold on

%Plot plane of B = 0
lengt = length(x1);
CO(:,:,1) = ones(lengt); % red
CO(:,:,2) = zeros(lengt); % green
CO(:,:,3) = zeros(lengt); % blue
mesh(zeros(lengt)+11.57,x2,x1/8,CO)
mesh(zeros(lengt)-11.57,x2,x1/8,CO)
CO(:,:,1) = zeros(lengt); % red
CO(:,:,2) = ones(lengt); % green
CO(:,:,3) = zeros(lengt); % blue
mesh(x1,x2,zeros(lengt),CO)

xlim([-15 15])
ylim([-30 30])
zlim([-2 5])
xlabel('Flapwise Tip Displacement')
ylabel('Flapwise Tip Velocity')
title('Flapwise Blade Bending Barrier Function, Green = operating region, Red = unsafe region')

  % Contour Plot
%plot unsafe and operation regions
x_u_1 = nsidedpoly(4, 'Center', [11.57+30 0], 'Sidelength', 60);
x_u_2 = nsidedpoly(4, 'Center', [-11.57-30 0], 'Sidelength', 60);
%original = 
g_x0 = @(x1,x2) 0.009531*x1.^2 + 0.001985*x1.*x2 + 0.00174*x2.^2 - 1;

g_x0 = @(x1,x2) -1 + 0.009541*x1.^2 + 0.002001*x1.*x2 + 0.001726*x2.^2;
 
 
% % BAD (no offset) on unsymmetric case ellipse = g_x0 = @(x1,x2) -1 + 0.01737*x1^2 + 0.003619*x1*x2 + 0.003171*x2^2;
% % 'GOOD' (with offset) on unsymmetric case ellispse = 
% g_x0 = @(x1,x2) -1 + 0.004816*x2.*(x1 - 1.5) + 0.02536*(x1 - 1.5).^2 + 0.004286*x2.^2;
% 
% % Union of inequalities example
% g_x0 = @(x1,x2) -1 + 0.008964*x1.^2 + 0.004039*x1.*x2 + 0.000561*x2.^2;
% -1 + 0.009531*x1^2 + 0.001986*x1*x2 + 0.00174*x2^2 - instead?????
 
figure
%Plot Solution
xx1 = -30:2:30;
xx2 = -30:2:30;
[x1,x2] = meshgrid(xx1,xx2);

% %Original
% B = 0.009531*x1.^2 + 0.001985*x1.*x2 + 0.00174*x2.^2 - 1.1;
% 
% % Offset (-1.5) barrier 
% B = 0.02054*x1.^2 + 0.003901*x1.*x2 - 2.05e-5*x1 + 0.004243*x2.^2 - 3.281e-6*x2 - 1.343;
% 
% % Union of inequalities example
% B = 0.008964*x1.^2 + 0.004039*x1.*x2 + 1.439e-6*x1 + 0.000561*x2.^2 + 3.6e-7*x2 - 1.1;
%  
 
[U,V2] = gradient(-B(x1,x2)); 
%contour(x1,x2,B)
hold on
quiver(x1,x2,U,V2)
hold on
plot(x_u_1, 'FaceColor', 'r')
hold on
plot(x_u_2,'FaceColor', 'r')
hold on
 
fimplicit(g_x0,'MeshDensity',300,'Color','g')
% [-15 15 -25 25]
xlim([-15 15])
ylim([-30 30])
xlabel('Flapwise Tip Displacement (m)')
ylabel('Flapwise Tip Velocity (m/s)')
title('Flapwise Blade Bending State Space, Green = operating region, Red = unsafe region')
