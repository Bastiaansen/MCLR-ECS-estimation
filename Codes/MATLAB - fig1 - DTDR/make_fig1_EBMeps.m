% This script uses the energy balance model with ocean heat uptake
% (EBM-epsilon) to find a good fit to DT and DR and uses that to compute
% the equilibrium warming. Implementation of the method detailed in the
% paper [Geoffroy et al, 2013].

%% Explicit times are needed:
t = linspace(1.5, length(DT)+0.5, length(DT))';
t_plot = linspace(0,10000,10000);

%% Part I: initial values (i.e. using epsilon = 1)

% Initial guessing (linear gregory fit)
coeff = polyfit(DT(1:150),DR(1:150),1);
F = coeff(2);
lambda = - coeff(1);
T_eq = F / lambda;

% Parameter computations
coeff = polyfit(t(30:150), log(1 - DT(30:150) / T_eq),1);
a_s = exp(coeff(2));
tau_s = - 1 / coeff(1);
a_f = 1 - a_s;


for i = 1:10
    tau_fs(i) = t(i) / (log(a_f) - log(1 - DT(i)/T_eq - a_s * exp(-t(i)/tau_s)));
end

tau_f = mean(tau_fs);

C = lambda / (a_f / tau_f + a_s / tau_s);
C0 = lambda * (tau_f * a_f + tau_s * a_s) - C;
gamma = C0 / (tau_f * a_s + tau_s * a_f);
phi_f = 1/a_f * ( 1 / (1 - tau_s/tau_f));
phi_s = (1 - a_f * phi_f)/a_s;

%% Computation of estimated temperature time series

Ts = T_eq - a_f * T_eq * exp(-t / tau_f) - a_s * T_eq * exp(-t / tau_s);
T0s = T_eq - a_f * phi_f * T_eq * exp(-t / tau_f) - a_s * phi_s * T_eq * exp(-t / tau_s);

%% Part II: the epsilon part in the EBM-epsilon thingy



epsilon = 1;
% Compute ocean uptake time series
Hs = gamma * (Ts - T0s);

%% Some iterations are needed to find the correct fitted values
for iter = 1:10

      

    % Regression:    
    i = i_end;

    Y = [DR(1:i)];
    [n,d] = size(Y);
    X = [ones(n,1), DT(1:i), Hs(1:i)];
    BETA = mvregress(X,Y);
    F = BETA(1);
    lambda = - BETA(2);
    epsilon = 1 - BETA(3);
    
    T_eq = F / lambda;
    T_eqs(i) = T_eq;
    
    % Parameter computations

    coeff = polyfit(t(30:150), log(1 - DT(30:150) / T_eq),1);
    a_s = exp(coeff(2));
    tau_s = - 1 / coeff(1);
    a_f = 1 - a_s;


    for i = 1:10
        tau_fs(i) = t(i) / (log(a_f) - log(1 - DT(i)/T_eq - a_s * exp(-t(i)/tau_s)));
    end

    tau_f = mean(tau_fs);

    C = lambda / (a_f / tau_f + a_s / tau_s);
    C0 = lambda * (tau_f * a_f + tau_s * a_s) - C;
    gamma = C0 / (tau_f * a_s + tau_s * a_f);
    phi_f = 1/a_f * ( 1 / (1 - tau_s/tau_f));
    phi_s = (1 - a_f * phi_f)/a_s;
    
    Ts = T_eq - a_f * T_eq * exp(-t / tau_f) - a_s * T_eq * exp(-t / tau_s);
    T0s = T_eq - a_f * phi_f * T_eq * exp(-t / tau_f) - a_s * phi_s * T_eq * exp(-t / tau_s);
    
    
    
    gamma = gamma / epsilon; % (Since above, we essentially have gamma'...)
    
    Hs = gamma * (Ts - T0s);
end
    
%% Compute values for plot
Ts = T_eq - a_f * T_eq * exp(-t_plot / tau_f) - a_s * T_eq * exp(-t_plot / tau_s);
T0s = T_eq - a_f * phi_f * T_eq * exp(-t_plot / tau_f) - a_s * phi_s * T_eq * exp(-t_plot / tau_s);
Hs = gamma * (Ts - T0s);

Rs = F - lambda * Ts - (epsilon - 1) * Hs;

%% Make plot

blueViolet = [54.1, 16.9, 88.6]/100;

plot(Ts, Rs, 'color', blueViolet, 'linewidth', 2.0)
scatter(T_eq, 0, '*',  'MarkerEdgeColor', blueViolet)
