%% Start with clean slate

close all
clear all

%% Paths

Path = '';
control_name = '_MPIESM12_control_1237.nc';
experiment_name = '_MPIESM12_abrupt4x_1000.nc';


%% Load-in data

ctrl_T = ncread([Path 'tas' control_name], 'tas');
ctrl_rlut = ncread([Path 'rlut' control_name], 'rlut');
ctrl_rsut = ncread([Path 'rsut' control_name], 'rsut');
ctrl_rsdt = ncread([Path 'rsdt' control_name], 'rsdt');

% Compute imbalance
ctrl_imb = ctrl_rsdt - ctrl_rsut - ctrl_rlut;

abr_T = ncread([Path 'tas' experiment_name], 'tas');
abr_rlut = ncread([Path 'rlut' experiment_name], 'rlut');
abr_rsut = ncread([Path 'rsut' experiment_name], 'rsut');
abr_rsdt = ncread([Path 'rsdt' experiment_name], 'rsdt');

% Compute imbalance
abr_imb = abr_rsdt - abr_rsut - abr_rlut;

%% Compute anomalies

% Initial values as mean of control experiment
T0 = mean(ctrl_T);
R0 = mean(ctrl_imb);

% anomalies as difference w.r.t. initial values
DT = abr_T - T0;
DR = abr_imb - R0;

DT = (DT(1:end-1)+DT(2:end))/2;
DR = (DR(1:end-1)+DR(2:end))/2;

%% Make figure

h = figure('Units','normalized','Position', [0.1 0 0.6 0.6]);

scatter(DT,DR, 'ro')
hold on

xlabel('$\Delta T$', 'Interpreter', 'latex')
ylabel('$\Delta R$', 'Interpreter', 'latex')


%% DOUBLE GREGORY
f = @(x,xdata) (x(3)+x(1)*xdata) .* (xdata < - (x(4)-x(3))/(x(2)-x(1))) + ...
    (x(4)+x(2)*xdata) .* (xdata >= - (x(4)-x(3))/(x(2)-x(1)));

coeff = polyfit(DT(1:10),DR(1:10),1);
lambda_1 = coeff(1);
f_1 = coeff(2);

x = lsqcurvefit(f,[lambda_1, lambda_1 + 0.1, f_1, f_1], DT, DR);
ECS_double_greg = - x(4)/x(2);

plot(DT, f(x, DT), 'g-')
scatter(ECS_double_greg, 0, 'g*')

%% TRIPLE GREGORY
f3 = @(x,xdata) (x(4)+x(1)*xdata) .* (xdata < - (x(5)-x(4))/(x(2)-x(1))) + ...
    (x(5)+x(2)*xdata) .* (xdata >= - (x(5)-x(4))/(x(2)-x(1))) .* (xdata < - (x(6)-x(5))/(x(3)-x(2))) + ...
    (x(6)+x(3)*xdata) .* (xdata >= - ((x(6)-x(5))/(x(3)-x(2))));

coeff = polyfit(DT(1:10),DR(1:10),1);
lambda_1 = coeff(1);
f_1 = coeff(2);

x = lsqcurvefit(f3,[lambda_1, lambda_1 + 0.5, lambda_1+1, f_1,f_1-1, f_1-3], DT, DR);
ECS_triple_greg = - x(6)/x(3);

plot(DT, f3(x, DT), 'b-')
scatter(ECS_triple_greg, 0, 'b*')
