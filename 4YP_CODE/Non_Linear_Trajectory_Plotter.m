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
% feta_beta - Blade-pitch [deg] (actually 90 minus pitch)
% omega_beta - Blade-pitch rate [deg/s] (actually minus pitch rate)

lambda_r = 2.025; % Max rotor velocity [rad/s]
lambda_delta = 441.42e-3; % Ultimate load limit of drive train torsion [rad]
%rated = 1.25!

dt = 0.01;
tspan = 0:dt:30;
x0 = [1.5 1.5*97 0 90 0];

[t,y] = ode45(@sys,tspan,x0);

 leny = length(y);
% plot(y(1:leny,1),y(1:leny,2))

figure
plot(t,y(1:leny,1))
hold on
%plot(t,zeros(1,leny)+441.42e-3)
xlabel('Time (s)')
ylabel('Torsion Angle (rad)')
title('Torsion Angle vs time during shutdown')

%% Pitch vs TSR plot
R = 61.5;
v_wind = 10;
Cq = zeros(1,leny);
Cq_pk = 0.08;
Cq_intcpt = 60;

for i = 1:leny
    if y(i,1)>0 && y(i,1)<10*v_wind/R

        TSR_i = R*y(i,1)/v_wind;
        Pitch_i = 90 - y(i,4);
        Cq(i) = (-((TSR_i/5)-1).^2 + 1).*(1-(Pitch_i/Cq_intcpt))*Cq_pk;
        else
        Cq(i) = 0;
    end
end


TSR = (0:10);
pitch = 0:10:100;

[TSR,pitch] = meshgrid(TSR,pitch);
c = (-((TSR/5)-1).^2 + 1).*(1-(pitch/Cq_intcpt))*Cq_pk;
figure
mesh(TSR,pitch,c)
hold on
[TSR_2,pitch_2] = meshgrid(10:20,0:10:100);
mesh(TSR_2,pitch_2,zeros(length(TSR_2)))
hold on
% [TSR_3,pitch_3] = meshgrid(-10:0,0:10:100);
% mesh(TSR_3,pitch_3,zeros(11))
% hold on
%plot3(y(1:leny,1)*R/v_wind,90-y(1:leny,4),Cq+0.002,'r')
hold on
%plot3(y(1,1)*R/v_wind,90-y(1,4),Cq(1)+0.002,'mo')
hold on
%plot3(y(leny,1)*R/v_wind,90-y(leny,4),Cq(leny)+0.002,'rx')
hold on
%Plot boundary for TSR (based on omega r)
lengt = length(pitch);
CO(:,:,1) = ones(lengt); % red
CO(:,:,2) = zeros(lengt); % green
CO(:,:,3) = zeros(lengt); % blue
mesh(zeros(lengt)+2.025*R/v_wind,pitch,TSR/100-0.05,CO)

xlabel('TSR')
ylabel('Pitch (degrees)')
zlabel('Torque Coefficient (approx)')
title('Trajectory of system across mesh of Torque Coefficient')
xlim([0 20])
ylim([0 100])



