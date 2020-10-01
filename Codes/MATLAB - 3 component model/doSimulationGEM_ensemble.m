%% doSimulationGEM_ensemble
% This script performs an ensemble of the three component global energy
% budget model in the script 'ThreeComponentGEM.m'. Subsequently, the
% obtained outputted time series on temperature T, albedo ALB and
% emissivity EMM are used to estimate the real model equilibrium warming
% and these estimates are compared to the real value to assess the power of
% different estimation techniques.
%
% MODEL:
% C_T T' = Q0 (1 - alpha) - sigma T^4 + mu + nu * W'
% alpha' = - eps_A * (alpha - alpha_0(T))
% sigma' = - eps_E * (simga - sigma_0(T))

%% Start with clean slate

clear all
close all

%% Ensemble parameters
N = 100; % Number of simulations in ensemble

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
nu = 0;

%% Simulation options
EndTime = 2000;
tspan = linspace(0,EndTime, 1000);

%% Initial conditions
% Temperature is specified and corresponding albedo and emissivity values
% are computed

T0 = 289;
ALB0 = alpha_0(T0);
EMM0 = sigma_0(T0);

y0 = [T0;ALB0;EMM0];

%% CO2 FORCING

% Initial value for mu is computed from initial conditions, since
% necessarily dT/dt = 0 at the start
mu0 = EMM0 * T0.^4 - Q0 * (1 - ALB0);

% Forcing
A0 = 5.35;
mu = @(t) A0 * log(4 + 0.*t) + mu0; % Instantaneous Quadruppling


%% Simulations -- Ensemble runs

ensemble_rem_rel_err = cell(N,1);


for sim = 1:N
    
    display( [...
        'Now performing simulation number ' num2str(sim) ...
        ' out of ' num2str(N)])
    
    %% Simulation
    [t,T,ALB,EMM] = ThreeComponentGEM(C_T,Q0,sigma_0,alpha_0, mu,eps_A, eps_E, nu, tspan, y0);
    
    %% Data processing
    DT = T - T0;
    DALB = ALB - ALB0;
    DEMM = (EMM - EMM0)/sigma;
    DR = Q0 .* (1 - ALB) + mu(t) - EMM .* T.^4;
    
    % Derivatives taken as forward differences
    DALBd = diff(DALB)./diff(t);
    DEMMd = diff(DEMM)./diff(t);
    
    % Use mean of subsequent points (so have same lengths as derivatives)
    DT = (DT(1:end-1)+DT(2:end))/2;
    DR = (DR(1:end-1)+DR(2:end))/2;
    DALB=(DALB(1:end-1)+DALB(2:end))/2;
    DEMM=(DEMM(1:end-1)+DEMM(2:end))/2;
    
    
    %% Compute the real equilibrium
    opts1=  optimset('display','off');
    DT_eq_real = fsolve(@(x) Q0*(1-alpha_0(x))+mu(t(end))-sigma_0(x).*x.^4,  T0, opts1) - T0;

    %% Perform the estimation techniques
    [estimates, estimates_info, ~] = perfom_estimations(DT,DR,DALB,DEMM,DALBd,DEMMd,C_T);
    
    %% Computation of (remaining) relative errors
    
    errors_rel = cell(length(estimates),1);
    errors_rel_rem = cell(length(estimates),1);
    
    for i = 1:length(estimates)
        estim = estimates{i};
        error_rel = abs( (estim - DT_eq_real) ./ estim );
        errors_rel{i} = error_rel;
        
        % Remaining error is maximum over later errors
        for j = 1:length(error_rel)
            error_rel_rem(j) = max(error_rel(j:end));
        end
        
        errors_rel_rem{i} = error_rel_rem;
    
    end
    
    ensemble_rem_rel_err{sim} = errors_rel_rem;    
    
end


%% Statistics
% Compute for each estimation technique the mean, standard deviation and 5
% and 95 percentile values

METHs = length(ensemble_rem_rel_err{1});

means = cell(METHs,1);
stds = cell(METHs,1);
perc5s = cell(METHs,1);
perc95s = cell(METHs,1);


% Make a loop over all METHODS (not over simulations/realisations!)
for i = 1:METHs
    
    M = length(ensemble_rem_rel_err{1}{i});
    
    mean_meth = NaN(1,M);
    std_meth = NaN(1,M);
    perc5_meth = NaN(1,M);
    perc95_meth= NaN(1,M);
    
    % Loop over the entries of one particular method
    for j = 1:M
        
       % Collect data from all simulations/realisations
       errs = NaN(1,N);
       for k = 1:N
          errs(k) = ensemble_rem_rel_err{k}{i}(j);           
       end
       
       % Perform statistics
       mean_meth(j) = mean(errs);
       std_meth(j) = std(errs);
       perc5_meth(j) = prctile(errs,5);
       perc95_meth(j)= prctile(errs,95);
        
    end
    
    means{i} = mean_meth;
    stds{i} = std_meth;
    perc5s{i} = perc5_meth;
    perc95s{i}= perc95_meth;
    
end




%% Plotting

figure('Units','normalized','Position', [0.1 0 0.6 0.6]);

colors = ['r', 'b', 'g', 'c', 'k', 'm'];

for i = 1:METHs
   loglog(t(2:end), means{i}, [colors(i) '-'], 'linewidth', 1.0)
   hold on
   loglog(t(2:end), perc5s{i}, [colors(i) '--'], 'linewidth', 1.0)
   loglog(t(2:end), perc95s{i},[colors(i) '--'], 'linewidth', 1.0)
end

xlabel('$t$', 'Interpreter', 'latex')
ylabel('$e_{rem}^{rel}$', 'Interpreter', 'latex')
axis([ 1 10^3 10^(-5) 10])































