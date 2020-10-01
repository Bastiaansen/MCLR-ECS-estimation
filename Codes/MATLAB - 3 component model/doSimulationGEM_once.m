%% doSimulationGEM_once
% This script performs one simulation of the three component global energy
% budget model in the script 'ThreeComponentGEM.m'. Subsequently, the
% obtained output time series on temperature T, albedo alb, and emissivity
% emm are used to estimate the real equilibrium warming and these estimates
% are compared to the real value.
%
% MODEL:
% C_T T' = Q0 (1 - alpha) - sigma T^4 + mu + nu * W'
% alpha' = - eps_A * (alpha - alpha_0(T))
% sigma' = - eps_E * (simga - sigma_0(T))

%% Start with clean slate
clear all
close all

%% Model parameters

% incoming solar radiation
Q0 = 341.3;

% Temperature dependent equilibrium values for albedo
alpha_1 = 0.7;
alpha_2 = 0.289;
T_1 = 260;
T_2 = 289;
K_A = 0.1;
alpha_0 = @(T) alpha_1 + (alpha_2 - alpha_1) * (1 + tanh(K_A * (T - (T_1+T_2)/2)))/2;
alpha_0d= @(T) (alpha_2-alpha_1)/2 * K_A * sech(K_A * (T - (T_1+T_2)/2)).^2;

% Temperature dependent equilibrium values for emissivity
sigma = 5.67 * 10^(-8);
sigma_1 = 0.7 * sigma;
sigma_2 = 0.6 * sigma;
K_E = 0.05;
Temm0 = 288;
sigma_0 = @(T) sigma_1 + (sigma_2-sigma_1) * (1 + tanh(K_E * (T - Temm0)))/2;
sigma_0d =@(T) (sigma_2-sigma_1)/2 * K_E * sech(K_E * (T-Temm0)).^2;

% Heat capacity
C_T = 10;

% Relaxation rates of albedo and emissivity
eps_A = 0.05;
eps_E = 0.1;

% Variance of noise
nu = 0.5;

%% Simulation options
EndTime = 200;
tspan = linspace(0,EndTime,1000);

%% Initial Conditions
% Temperature is specified and corresponding albedo and emissivity are
% computed

T0 = 289;
ALB0 = alpha_0(T0);
EMM0 = sigma_0(T0);

y0 = [T0; ALB0; EMM0];

%% CO2 forcing

% Initial value for mu is computed from initial conditions, since
% necessarily dT/dt = 0 at the start.
mu0 = EMM0 * (T0).^4 - Q0 * (1 - ALB0);

% Forcing
A0 = 5.35;
mu = @(t) A0 * log(4 + 0.*t) + mu0; % Instantaneous Quadruppling

%% Run the simulation
[t,T,ALB,EMM] = ThreeComponentGEM(C_T,Q0,sigma_0,alpha_0, mu,eps_A, eps_E, nu, tspan, y0);

%% Data processing
DT = T - T0;
DALB = ALB - ALB0;
DEMM = (EMM - EMM0)/sigma;
DR = Q0 .* (1 - ALB) + mu(t) - EMM .* T.^4;

% Derivaties taking as forward differences
DALBd = diff(DALB)./diff(t);
DEMMd = diff(DEMM)./diff(t);

% Use mean of subsequent points (so have same length as derivatives)
DT = (DT(1:end-1)+DT(2:end))/2;
DR = (DR(1:end-1)+DR(2:end))/2;
DALB=(DALB(1:end-1)+DALB(2:end))/2;
DEMM=(DEMM(1:end-1)+DEMM(2:end))/2;

%% Compute the real equilibrium
opts1=  optimset('display','off');
DT_eq_real = fsolve(@(x) Q0*(1-alpha_0(x))+mu(t(end))-sigma_0(x).*x.^4,  T0, opts1) - T0;


%% Perform the estimation techniques
[estimates, estimates_info, estimates_eigenvalues] = perfom_estimations(DT,DR,DALB,DEMM,DALBd,DEMMd,C_T);


%% Plotting of the equilibrium warming predictions

h_EQ = figure('Units', 'normalized', 'Position', [0.1 0 0.6 0.6]);
hold on
plot([t(1) t(end)], [DT_eq_real, DT_eq_real], 'k:', 'linewidth', 2.0)
plot(t(2:end), estimates{1}, 'r-', 'linewidth', 2.0)
plot(t(2:end), estimates{2}, 'b-', 'linewidth', 2.0)
plot(t(2:end), estimates{3}, 'g-', 'linewidth', 2.0)
plot(t(2:end), estimates{4}, 'c-', 'linewidth', 2.0)
plot(t(2:end), estimates{5}, 'k-', 'linewidth', 2.0)
plot(t(2:end), estimates{6}, 'm-', 'linewidth', 2.0)
xlabel('$t$', 'Interpreter', 'latex', 'fontsize', 50)
ylabel('$\Delta T_*$', 'Interpreter', 'latex', 'fontsize', 50)
title('Equilibrium warming', 'Interpreter', 'latex', 'fontsize', 50)
axis([t(1) t(end) 0 ceil(DT_eq_real)+1])
legend('True value', 'Raw Simulation Value', ...
    'Gregory', 'Double Gregory', 'System Fit [T,ALB]', ...
    'System Fit [T,EMM]', 'System Fit [T,ALB,EMM]')

%% Computation of errors
for i = 1:length(estimates)
    estim = estimates{i};
    error_rel = abs( (estim - DT_eq_real)./estim);
    errors_rel{i} = error_rel;
    
    % Remaining error is maximum over later errors
    for j=1:length(error_rel)
        error_rel_rem(j) = max(error_rel(j:end));
    end
    
    errors_rel_rem{i} = error_rel_rem;
    
end

%% Plotting of remaining relative errors

h_err = figure('Units', 'normalized', 'Position', [0.1 0 0.6 0.6]);
hold on
plot(t(2:end), errors_rel_rem{1}, 'r-', 'linewidth', 2.0)
plot(t(2:end), errors_rel_rem{2}, 'b-', 'linewidth', 2.0)
plot(t(2:end), errors_rel_rem{3}, 'g-', 'linewidth', 2.0)
plot(t(2:end), errors_rel_rem{4}, 'c-', 'linewidth', 2.0)
plot(t(2:end), errors_rel_rem{5}, 'k-', 'linewidth', 2.0)
plot(t(2:end), errors_rel_rem{6}, 'm-', 'linewidth', 2.0)
xlabel('$t$', 'Interpreter', 'latex', 'fontsize', 50)
ylabel('Remaining relative error', 'Interpreter', 'latex', 'fontsize', 50)
title('error in estimates', 'Interpreter', 'latex', 'fontsize', 50)
legend('True value', 'Raw Simulation Value', ...
    'Gregory', 'Double Gregory', 'System Fit [T,ALB]', ...
    'System Fit [T,EMM]', 'System Fit [T,ALB,EMM]')
axis([t(1) t(end) 0 1])

%% Computation of the real eigenvalues
% real equilibrium values, Jacobian and eigenvalues
T_star = T0 + DT_eq_real;
ALB_star = alpha_0(T_star);
EMM_star = sigma_0(T_star)/sigma;
J_real = [...
    -4 * sigma * EMM_star * T_star^3 / C_T, -Q0/C_T, -sigma * T_star^4 / C_T; ...
    eps_A * alpha_0d(T_star), -eps_A, 0; ...
    eps_E * sigma_0d(T_star)/sigma, 0, -eps_E];
lambda_real = eig(J_real);

%% Plotting of eigenvalues

h_EV = figure('Units', 'normalized', 'Position', [0.1 0 0.6 0.6]);
hold on
plot([t(1) t(end)], [lambda_real, lambda_real], 'k:', 'linewidth', 2.0)
plot(t(2:end), estimates_eigenvalues{2}, 'b-', 'linewidth', 2.0)
plot(t(2:end), estimates_eigenvalues{3}, 'g-', 'linewidth', 2.0)
plot(t(2:end), estimates_eigenvalues{4}, 'c-', 'linewidth', 2.0)
plot(t(2:end), estimates_eigenvalues{5}, 'k-', 'linewidth', 2.0)
plot(t(2:end), estimates_eigenvalues{6}, 'm-', 'linewidth', 2.0)
xlabel('$t$', 'Interpreter', 'latex')
ylabel('$\lambda$', 'Interpreter', 'latex')
legend('True', 'True', 'True', ...
    'Gregory', 'Double Gregory', 'Double Gregory', ...
    'System Fit [T,ALB]', 'System Fit [T,ALB]', ...
    'System Fit [T,EMM]', 'System Fit [T,EMM]', ...
    'System Fit [T,ALB,EMM]', 'System Fit [T,ALB,EMM]', 'System Fit [T,ALB,EMM]')
axis( [t(1) t(end) -0.5 0.1] )

