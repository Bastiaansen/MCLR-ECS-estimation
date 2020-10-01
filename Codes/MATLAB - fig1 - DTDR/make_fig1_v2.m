% This script makes 'figure 1': The Gregory plot DT vs DR that shows
% nonlinearity of the plot, and visualizes the various estimation
% techniques including their estimated warmings.

% Paths below need to be directed to the data. For our plots, we have used
% LongRunMIP data for CESM 1.0.4. (N.B. not tested excessively on other data)

%% Start with clean slate
close all
clear all

%% Paths

Path = 'CESMdata/';
control_name = '_CESM104_control_1000.nc';
experiment_name = '_CESM104_abrupt4x_5900.nc';


%% Load-in data + computations

ctrl_T = ncread([Path 'tas' control_name], 'tas');
ctrl_rlut = ncread([Path 'rlut' control_name], 'rlut');
ctrl_rsut = ncread([Path 'rsut' control_name], 'rsut');
ctrl_rsdt = ncread([Path 'rsdt' control_name], 'rsdt');

% Compute imbalance, albedo & emissivity
ctrl_imb = ctrl_rsdt - ctrl_rsut - ctrl_rlut;
ctrl_alb = ctrl_rsut ./ ctrl_rsdt;
ctrl_emm = ctrl_rlut ./ (ctrl_T.^4);

abr_T = ncread([Path 'tas' experiment_name], 'tas');
abr_rlut = ncread([Path 'rlut' experiment_name], 'rlut');
abr_rsut = ncread([Path 'rsut' experiment_name], 'rsut');
abr_rsdt = ncread([Path 'rsdt' experiment_name], 'rsdt');

% Compute imbalance, albedo & emissivity
abr_imb = abr_rsdt - abr_rsut - abr_rlut;
abr_alb = abr_rsut ./ abr_rsdt;
abr_emm = abr_rlut ./ (abr_T.^4);

%% Compute anomalies

% Initial values as means of control experiment
T0 = mean(ctrl_T);
R0 = mean(ctrl_imb);
ALB0 = mean(ctrl_alb);
EMM0 = mean(ctrl_emm);

% Anomalies as difference w.r.t. initial values
DT = abr_T - T0;
DR = abr_imb - R0;
DALB = abr_alb - ALB0;
DEMM = abr_emm - EMM0;

%% Compute derivatives (and compute in-between values)

ALBd = diff(DALB);
EMMd = diff(DEMM);

DT = (DT(1:end-1)+DT(2:end))/2;
DR = (DR(1:end-1)+DR(2:end))/2;
DALB = (DALB(1:end-1)+DALB(2:end))/2;
DEMM = (DEMM(1:end-1)+DEMM(2:end))/2;

%% Make figure -- scatter plot

h = figure('Units','normalized','Position', [0.1 0 0.6 0.6]);

% Create size of scatter elements based on amount of years
scat_size = [linspace(10,0.1,1000), 0.1*ones(1,length(DT)-1000)];
k = 20;
scat_size = 0.2 + 15 * exp(- k * [1:length(DT)]/length(DT));

scatter(DT,DR, scat_size, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'MarkerFaceAlpha', .5, 'MarkerEdgeAlpha', .2)
hold on

xlabel('$\Delta T$', 'Interpreter', 'latex')
ylabel('$\Delta R$', 'Interpreter', 'latex')

%% Mark the best estimated range (hard-coded here!)
DT_range = [6.73, 6.80];
best_range_plot = plot(DT_range, [0,0], 'color', [17 17 17]/255, 'linewidth', 5);
best_range_plot.Color(4) = 0.5;

%% Option for fitting
i_end = 300; % Last year to use in the fits

%% GREGORY METHODS
make_fig1_gregoryFits

%% 3-EXP fits [Proistosecu and Huybers, 2017]
make_fig1_3EXP

%% EBM-epsilon fits [Geoffroy et al, 2013]
make_fig1_EBMeps

%% System fit method
make_fig1_sysFit

%% Make-up of plot

axis([1 7 -1 6])
plot([0 8], [0 0], 'k:')

legend()