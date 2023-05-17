function dydt = sys_test(t,y)

%% SOS tools tester system
d = 1.2;
dydt = [ y(2)
    (-y(1) + (d*(y(1).^3)/3) -y(2))];





