%% SOS Lyapunov Tester! - Flapwise First


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


syms x1 x2;

% =============================================
% First, initialize the sum of squares program
prog = sosprogram([x1;x2]);

f = [x2; (1/M_flap)*(F_aero - K_flap*(x1) - B_flap*x2)];

% =============================================
% The Lyapunov function V(x): 
[prog,V] = sospolyvar(prog,[x1^2; x2^2;x1*x2],'wscoeff');

% =============================================
% Next, define SOSP constraints

%Declare a decision variable t (NEW)
%prog = sosdecvar(prog, t);

% Constraint 1 : V(x) - (x1^2 + x2^2 + x3^2) >= 0
prog = sosineq(prog,V-(x1^2+x2^2));

% Constraint 2: -dV/dx*(x3^2+1)*f >= 0
expr = -(diff(V,x1)*f(1)+diff(V,x2)*f(2));

prog = sosineq(prog,expr);

% Set Objective function? (new bit)
%prog = sossetobj(prog, -t); %maximise t? so that gradV is the most positive?

% =============================================
% And call solver
solver_opt.solver = 'sedumi';
prog = sossolve(prog,solver_opt);

% =============================================
% Finally, get solution
SOLV = sosgetsol(prog,V)

%%Plot Solution
xx1 = -30:0.5:30;
xx2 = -30:0.5:30;
[x1,x2] = meshgrid(xx1,xx2);

V = matlabFunction(SOLV-780);

 
figure
mesh(x1,x2,V(x1,x2)+780)
hold on

%Plot plane of B = 0
lengt = length(x1);
CO(:,:,1) = ones(lengt); % red
CO(:,:,2) = zeros(lengt); % green
CO(:,:,3) = zeros(lengt); % blue
mesh(zeros(lengt)+11.57,x2,x1*90,CO)
mesh(zeros(lengt)-11.57,x2,x1*90,CO)
CO(:,:,1) = zeros(lengt); % red
CO(:,:,2) = ones(lengt); % green
CO(:,:,3) = zeros(lengt); % blue
mesh(x1,x2,zeros(lengt)+780,CO)


xlim([-15 15])
ylim([-30 30])
zlim([-500 5000])
xlabel('Flapwise Tip Displacement (m)')
ylabel('Flapwise Tip Velocity (m/s)')
zlabel('Lyapunov Function')
%title('Flapwise Blade Bending Lyapunov Function')

%%Quiver plot
xx1 = -30:4:30;
xx2 = xx1;

[x1_grad,x2_grad] = meshgrid(xx1,xx2);

u = zeros(size(x1_grad));
v = zeros(size(x2_grad));

% we can use a single loop over each element to compute the derivatives at
% each point (y1, y2)
t=0; % we want the derivatives at each point at t=0, i.e. the starting time
f = @(t,y) [ y(2); (1/M_flap)*(F_aero - K_flap*(y(1)) - B_flap*y(2))];

for i = 1:numel(x1_grad)
    Yprime = f(t,[x1_grad(i); x2_grad(i)]);
    u(i) = Yprime(1);
    v(i) = Yprime(2);
end


%plot unsafe and operation regions
x_u_1 = nsidedpoly(4, 'Center', [11.57+30 0], 'Sidelength', 60);
x_u_2 = nsidedpoly(4, 'Center', [-11.57-30 0], 'Sidelength', 60);

figure
% 
% fimplicit(V,[-30 30 -30 30],'g')

contour(x1,x2,V(x1,x2)-1000,'g')
hold on
quiver(xx1,xx2,u,v,'b')
hold on

plot(x_u_1, 'FaceColor', 'r')
hold on
plot(x_u_2,'FaceColor', 'r')

xlim([-15 15])
ylim([-30 30])
xlabel('Flapwise Tip Displacement (m)')
ylabel('Flapwise Tip Velocity (m/s)')




