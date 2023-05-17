%% Plotting shutdown pitch model!
clear all
close all
% 
% f = [y(2);
%     -0.6*y(2) - 0.0894*y(1)];
p_0 = 10; % pitch at which 
ratio = 8/(90-p_0); %Ratio of spring term to damper term, 
%(to ensure that pitch acceleration is zero at operating pitch!)

zeta = 0.85;

lambda = -(2*zeta)^2*ratio;
%lambda = -0.1782;

A = [0 1;lambda*[ratio 1]];
A = [0 1;-0.0894 -0.6];
B = [0;0];
C = eye(2);
D = [0;0];
%% plug into initial
sys = ss(A,B,C,D);

% Run 'initial'
pitch_initial  = 10;

for pitch_initial = 0
str = 'Intial Pitch = ' + string(pitch_initial) + ' (degrees)';

x0 = [90-pitch_initial,0]; % Initial Condition

[y,t,x] = initial(sys,x0);

leny = length(y);

figure
plot(t(1:leny),90-y(1:leny,1))
hold on

plot([0 (90-pitch_initial)/8 20],[pitch_initial 90 90])
legend('Modelled Response','Ideal Response')
xlabel('Time (s)')
ylabel('Pitch (degrees)')

xlim([0 20])
ylim([0 110])
%title(str)
end
