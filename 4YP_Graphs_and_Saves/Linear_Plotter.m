%Program to look at bode plots etc... from linearised model
clear ; close all;
%load linearised output file 
[data] = ReadFASTLinear('IEA-15-240-RWT-Monopile.1.lin');

A = data.A;
B = data.B;
C = data.C;
D = data.D;

% Want - pitch angle(9), gen torq(8) to gen speed(10) and nacelle foreaft(19)
% For other format!!!! - pitch = 6, torq = 5, gen speed = 6,foreaft (m)= 7 
% Plug matrices into state space

sys = ss(A,B,C,D);
Func = tf(sys);

% Find individual transfer functions for each in/out
Ang_to_Speed = tf(Func.Numerator(10,9),Func.Denominator(10,9));

Ang_to_Foreaft = tf(Func.Numerator(19,9),Func.Denominator(19,9));

Torq_to_Speed = tf(Func.Numerator(10,8),Func.Denominator(10,8));

Torq_to_Foreaft = tf(Func.Numerator(19,8),Func.Denominator(19,8));

% Find individual transfer functions for each in/out
% Ang_to_Speed = ss(A,B(:,6),C(6,:),D(6,6));
% 
% Ang_to_Foreaft = ss(A,B(:,6),C(7,:),D(7,6));
% 
% Torq_to_Speed = ss(A,B(:,5),C(6,:),D(6,5));
% 
% Torq_to_Foreaft = ss(A,B(:,5),C(7,:),D(7,5));

%Plot Bode Plots of all of the above 
%Start with torq to speed and torq to foreaft

bode(Torq_to_Speed)
title('Generator Torque to Generator Speed')

% figure
% bode(Torq_to_Foreaft)
% title('Generator Torque to Fore-aft Acceleration')

figure
bode(Ang_to_Speed)
title('Blade Pitch to Generator Speed')
% 
% figure
% bode(Ang_to_Foreaft)
% title('Blade Pitch to Fore-aft Acceleration')

figure
pzmap(Torq_to_Speed)
title('Generator Torque to Generator Speed')

figure
pzmap(Ang_to_Speed)
title('Blade Pitch to Generator Speed')
eig(sys.A)

% atos_bal = balreal(sys_atos,1e-3) 
% pzmap(balreal(sys_atos,1e-4))
% pzmap(balreal(sys_atos,1e-3))
% size(balreal(sys_atos,1e-3))

% Find emptyinesss? of A and B 
% spy(sys.B)
% spy(sys.A)
% % Find particular states of interest before finding all transfer functions
%sys_atos = ss(A,B(9,:),C(:,10),D(9,10));
% Find eigen values eig(sys.A)

% Mod red - reduces system - removing the unwanted states from balreal?
% Bal real alligns gramians to kinda rank the the most observable and
% controllable states :))
% 6th output = azimuth

% figure
% bode(Ang_to_Speed)

% hold on; 
% bode(modred(balreal(Ang_to_Speed),23:40)) %Keep first 20 states only!
% 
% 
