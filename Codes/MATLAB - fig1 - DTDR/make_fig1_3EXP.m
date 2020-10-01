% This script uses the 3-EXP method to obtain fits for DT and DR to a
% function that is a sum of three exponentials. This method originates in
% the paper [Proistosecu and Huybers, 2017]

%% Define function & set options to disable additional information
opts1=  optimset('display','off');

f = @(x,t) ...
    [ x(6) * (1 - x(4)*exp(-t/x(1)) - x(5) * exp(-t/x(2)) - (1-x(4)-x(5)) * exp(-t/x(3))); ...
    x(9) * (x(7) * exp(-t/x(1)) + x(8) * exp(-t/x(2)) + (1-x(7)-x(8)) * exp(-t/x(3)))];
% Initial guess
x_init = [0.5, 25, 200, 0.33, 0.33, 4, 0.33, 0.33, 7];

%% FIT ON YEARS 1-i_end

x = lsqcurvefit(f, x_init, [1:1:length(DT(1:i_end))], [DT(1:i_end)';DR(1:i_end)'],[],[], opts1);
ECS = x(6);
F0 = x(9);

%% Make plot
ts = linspace(0,10000,1000);
y_res = f(x, ts);
plot(y_res(1,:),y_res(2,:), 'g-', 'linewidth', 2.0)
scatter(ECS, 0, 'g*')