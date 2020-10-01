%% This script makes 'figure 1': the Gregory plot DT vs DR that shown nonlinearity
% Paths below need to be directed to the data. For our plots, we have used
% LongRunMIP data for CESM1.0.4 model. (N.B. Not tested on other data!)

%% Start with clean slate

close all
clear all

%% Paths

Path = 'CESMdata/';
control_name = '_CESM104_control_1000.nc';
experiment_name = '_CESM104_abrupt4x_5900.nc';


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

%% Make figure

h = figure('Units','normalized','Position', [0.1 0 0.6 0.6]);

scatter(DT,DR, 'ro')
hold on

xlabel('$\Delta T$', 'Interpreter', 'latex')
ylabel('$\Delta R$', 'Interpreter', 'latex')

%% Linear Regression 1: whole dataset

coeff = polyfit(DT,DR,1);
x_intercept = coeff(2) / (-coeff(1));
y_intercept = coeff(2);

plot([0 2*x_intercept], [1,-1] * y_intercept, 'b-')

scatter(x_intercept, 0, 'b*')

%% Linear Regression 2: y 1-20

coeff = polyfit(DT(1:20),DR(1:20),1);
x_intercept = coeff(2) / (-coeff(1));
y_intercept = coeff(2);

plot([0 2*x_intercept], [1,-1] * y_intercept, 'k-')

% scatter(x_intercept, 0, 'k*')

%% Linear Regression 3: y 20-150

coeff = polyfit(DT(20:150),DR(20:150),1);
x_intercept = coeff(2) / (-coeff(1));
y_intercept = coeff(2);

plot([0 2*x_intercept], [1,-1] * y_intercept, 'g-')

scatter(x_intercept, 0, 'g*')



%% Linear Regression 5: y last 100 years

coeff = polyfit(DT(end-100:end),DR(end-100:end),1);
x_intercept = coeff(2) / (-coeff(1));
y_intercept = coeff(2);

plot([0 2*x_intercept], [1,-1] * y_intercept, 'y-')

scatter(x_intercept, 0, 'y*')


%% Make-up of plot

axis([0 8 -1 8])
plot([0 8], [0 0], 'k:')



