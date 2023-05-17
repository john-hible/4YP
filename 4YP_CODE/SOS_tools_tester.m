%% SOS Tools Tester

clear all
close all
%% MAKE SURE SOSTOOLs & SeDuMi (AND YALMIP?) ARE IN PATH

options.solver='SeDuMi';
% 
%% Formulate p as SOS
% syms x y;
% p = 2*x^4 + 2*x^3*y - x^2*y^2 + 5*y^4;
% [Q, Z] = findsos(p, options)
% % Returns matrices Q and Z so that Z.'*Q*Z = p

%% Optimisation problem
% 
% syms g p1 p2 p3; % defines variables to use!
% d = 4; % Degree bound for SOS hierarchy? higher is better but more to compute
% 
% obj = g; % Specify objective to MINIMISE
% 
% P = [p1 p2;p2 p3];
% 
% A = [-1 1;-2 -1];
% G = [1;1];
% C = [1 1];
% 
% M = [-A'*P-P*A-C'*C, -P*G; -G'*P, g];
% 
% 
% ineqs = [M,P]; % Specify all inequalities (>= 0)
% 
% [t,vars,xopt] = findbound(obj, ineqs, 2*d, options) % formulate in function
% % Returns t ~ 1.3911

%% Lyapunov Problem
% 
% pvar x1 x2;
% 
% %Construct the vector field (x' = f(x))
% f = [(-2*x1);(x1-x2)];
% 
% %Create SOS program
% prog = sosprogram([x1;x2]);
% 
% %Create Lyapunov Function V
% [prog, V] = sospolyvar(prog,[x1; x2; x1^2; x1*x2 ;x2^2]);
% 
% %Constrain so V is positive and gradV is negative 
% prog = sosineq(prog,V - (x1^2 + x2^2));
% 
% prog = sosineq(prog, - (diff(V,x1)*f(1) + diff(V,x2)*f(2)));
% 
% %Solve!! (solver already called at top)
% prog = sossolve(prog,options);
% 
% %Get Solution!
% SOLV = sosgetsol(prog,V)
% 
% %% Plot Resulting Lyapunov Function
% 
% x1 = -20:0.4:20;
% x2 = x1;
% x = [x1 ; x2];
% [xx1,xx2] = meshgrid(x1,x2);
% coeff = full(SOLV.coefficient);
% % V = zeros(length(x1),length(x2));
% 
% V = coeff(1)*xx1.^2 + coeff(2)*xx1.*xx2 + coeff(3)*xx2.^2; 
% % for i = 1:length(x1)
% %     for j = 1:length(x2)
% %         x = [x1(i); x2(j)];
% %         V(i,j) = x'*P*x;        
% %     end
% % end
% 
% [U,V2] = gradient(-V,0.4);
% contour(x1,x2,V)
% 
% hold on
% quiver(x1,x2,U,V2)
% hold off
% 
% figure
% mesh(x1,x2,V)

%% Barrier Function Problem!

syms x1 x2 d;

%Create SOS program
prog = sosprogram([x1;x2;d]);

%Construct the vector field (x' = f(x))

%d = 1;

f = [x2; -x1+d*(x1^3)/3-x2];

%Define initial ans unsafe spaces
r_x0 = 0.5;
g_x0 = ((x1-1.5)^2+x2^2-r_x0^2); % <0
%g_x0 = ((x1-1.5)^2+x2^2-0.25); % <0
g_xu = ((x1+1)^2+(x2+1)^2-0.16); % <0

g_xD = (d - 0.8)*(d - 1.2); %<0


% Define set of monomials for input to sos poly variables - INCLUDE 0th
% Order!!!
VEC_4 = monomials([x1;x2;d],0:4); % B up to order 4!
VEC_2 = monomials([x1;x2;d],0:3); %Sigmas up to order 2

%Create Barrier Function B, and sigma_1 and sigma_2
[prog, B] = sospolyvar(prog,VEC_4);

[prog,sig_1] = sossosvar(prog,VEC_2);
[prog,sig_2] = sossosvar(prog,VEC_2);
[prog,sig_3] = sossosvar(prog,VEC_2);
[prog,sig_4] = sossosvar(prog,VEC_2);
[prog,sig_5] = sossosvar(prog,VEC_2);

%Constrain so B is positive/negative in respective regions and gradB is negative 
prog = sosineq(prog,-B -0.1 + sig_1*g_x0 + sig_3*g_xD);
prog = sosineq(prog, B -0.1 + sig_2*g_xu + sig_4*g_xD);

prog = sosineq(prog, - (diff(B,x1)*f(1) + diff(B,x2)*f(2)) + sig_5*g_xD);

%
%% Lyapunov constraint!
% prog = sosineq(prog,B-(x1^2+x2^2));
% 
% prog = sosineq(prog, - (diff(B,x1)*f(1) + diff(B,x2)*f(2)));


%Solve!! (solver already called at top)
prog = sossolve(prog,options);

%Get Solution!
SOLV = sosgetsol(prog,B)
%return
xx1 = -5:0.2:5;
xx2 = xx1;
[x1,x2] = meshgrid(xx1,xx2);

%B =  -0.3361*x1.^4 - 1.0429*x1.^3.*x2 - 0.48887*x1.^2*x2.^2 - 0.42805*x1.*x2.^3 - 3.9818e-08*x2.^4 + 1.9479e-05*x1.^3 - 2.6885e-05*x1.^2.*x2 - 1.9126*x1.*x2.^2 - 0.6033*x2.^3 + 2.0166*x1.^2 + 3.1288*x1.*x2 + 4.2517*x2.^2 + 1.7244e-05*x1 + 6.0251e-06*x2 - 4.332;

%B = - 0.3361*x1.^4 - 1.043*x1.^3.*x2 + 1.944e-5*x1.^3 - 0.4889*x1.^2.*x2.^2 - 2.683e-5*x1.^2.*x2 + 2.017*x1.^2 - 0.4281.*x1.*x2.^3 - 1.913.*x1.*x2.^2 + 3.129.*x1.*x2 + 1.721e-5.*x1 - 3.918e-8.*x2.^4 - 0.6033.*x2.^3 + 4.252.*x2.^2 + 6.016e-6.*x2 - 4.332;
B = matlabFunction(SOLV-10);

%% Simulate and plot trajectories!
dt = 0.01;
tspan = 0:dt:30;
x0 = [1.9 0.29];

[t,y_1] = ode45(@sys_test,tspan,x0);

x0 = [1.2 -0.41];
[t,y_2] = ode45(@sys_test,tspan,x0);



%% Plot stuff
figure
mesh(x1,x2,B(0.8,x1,x2))
hold on

%Plot plane of B = 0
lengt = length(x1);
CO(:,:,1) = ones(lengt); % red
CO(:,:,2) = zeros(lengt); % green
CO(:,:,3) = zeros(lengt); % blue
mesh(x1,x2,zeros(lengt)-10,CO)

%plot unsafe and initial regions
x_0 = nsidedpoly(1000, 'Center', [1.5 0], 'Radius', r_x0);
x_u = nsidedpoly(1000, 'Center', [-1 -1], 'Radius', 0.4);

plot(x_0, 'FaceColor', 'g')
hold on
plot(x_u,'FaceColor', 'r')

xlabel('x1')
ylabel('x2')
title('Plot showing Barrier Function (red = <0, blue = >0) Green Circle = operating region, Purple Circle = unsafe region')

figure 

B_cont = B(0.8,x1,x2) + 10;

xx1 = -5:0.5:5;
xx2 = xx1;

[x1_grad,x2_grad] = meshgrid(xx1,xx2);

u = zeros(size(x1_grad));
v = zeros(size(x2_grad));

% we can use a single loop over each element to compute the derivatives at
% each point (y1, y2)
t=0; % we want the derivatives at each point at t=0, i.e. the starting time
d = 0.8;
f = @(t,y) [ y(2);(-y(1) + (d*(y(1).^3)/3) -y(2))];

for i = 1:numel(x1_grad)
    Yprime = f(t,[x1_grad(i); x2_grad(i)]);
    u(i) = Yprime(1);
    v(i) = Yprime(2);
end


contour(x1,x2,B_cont,'k--')
hold on
quiver(xx1,xx2,u,v,2.5)
hold on
x_0 = nsidedpoly(1000, 'Center', [1.5 0], 'Radius', r_x0);
x_u = nsidedpoly(1000, 'Center', [-1 -1], 'Radius', 0.4);
plot(x_0, 'FaceColor', 'g')
hold on
plot(x_u,'FaceColor', 'r')

hold on
plot(y_1(1:length(y_1),1),y_1(1:length(y_1),2),'b')

hold on
plot(y_2(1:length(y_2),1),y_2(1:length(y_2),2),'b')

xlim([-4 5])
ylim([-4 3])

xlabel('x_1')
ylabel('x_2')

