%% Repeat example from SOS guide but with maximise the area of operating region

clear all
close all
%% MAKE SURE SOSTOOLs & SeDuMi (AND YALMIP?) ARE IN PATH

options.solver='SeDuMi';

%% Barrier Function Problem!

syms x1 x2 q1 q2 q3 r;

%Create SOS program
prog = sosprogram([x1;x2]);
prog = sosdecvar(prog, [q1;q2;q3;r]);

%Define initial ans unsafe spaces
% g_x0 = q1*(x1)^2 + q2*x1*x2 + q3*(x2)^2 - 1;
g_x0 = q1*(x1-1.5)^2 + q2*x2^2 + r;
g_xu = ((x1+1)^2+(x2+1)^2-0.1); % <0

%Construct the vector field (x' = f(x))
f = [x2; -x1+(x1^3)/3-x2];

% Define set of monomials for input to sos poly variables - INCLUDE 0th
% Order!!!
VEC_4 = monomials([x1; x2],0:4); % B up to order 4!
VEC_2 = monomials([x1; x2],0:4); %Sigmas up to order 2

%Create Barrier Function B, and sigma_1 and sigma_2
[prog, B] = sospolyvar(prog,VEC_4);

[prog,sig_1] = sossosvar(prog,VEC_2);
[prog,sig_2] = sossosvar(prog,VEC_2);

%Constrain so B is positive/negative in respective regions and gradB is negative 
% Inequalities are >= 0 
prog = sosineq(prog,-B -0.1 + g_x0);
prog = sosineq(prog,r);
prog = sosineq(prog, B -0.1 + sig_2*g_xu);

prog = sosineq(prog, - (diff(B,x1)*f(1) + diff(B,x2)*f(2)));

%Define Objective Function
prog = sossetobj(prog,-r); %Minimise Tr(Q) = q1*q3

%Solve!! (solver already called at top)
prog = sossolve(prog,options);

%Get Solution!
SOLV = sosgetsol(prog,B)
R = (sosgetsol(prog,r))^0.5;
%return
xx1 = -20:0.5:20;
xx2 = xx1;
[x1,x2] = meshgrid(xx1,xx2);

B = - 1.757e-5*x1.^4 - 6.884e-6*x1.^3.*x2 + 1.658e-11*x1.^3 - 1.568e-6*x1.^2.*x2.^2 + 6.431e-12*x1.^2.*x2 + 0.0001054*x1.^2 - 1.806e-7*x1.*x2.^3 + 9.243e-10*x1.*x2.^2 + 2.066e-5*x1.*x2 - 3.127e-11*x1 - 9.733e-9*x2.^4 + 1.693e-10*x2.^3 + 6.71e-5*x2.^2 - 4.086e-12*x2 + 0.3975;
 
figure
mesh(x1,x2,B - 10)
hold on

%Plot plane of B = 0
lengt = length(x1);
CO(:,:,1) = ones(lengt); % red
CO(:,:,2) = zeros(lengt); % green
CO(:,:,3) = zeros(lengt); % blue
mesh(x1,x2,zeros(lengt)-10,CO)

%plot unsafe and initial regions
x_0 = nsidedpoly(1000, 'Center', [1.5 0], 'Radius', R);
x_u = nsidedpoly(1000, 'Center', [-1  -1], 'Radius', 0.1^0.5);
plot(x_0, 'FaceColor', 'g')
hold on
plot(x_u,'FaceColor', 'r')

xlabel('x1')
ylabel('x2')
title('Plot showing Barrier Function (red = <0, blue = >0) Green Circle = operating region, Purple Circle = unsafe region')

figure 

[U,V2] = gradient(-B,0.4);
contour(x1,x2,B)
hold on
% quiver(x1,x2,U,V2)
hold on
plot(x_0, 'FaceColor', 'g')
hold on
plot(x_u,'FaceColor', 'r')



