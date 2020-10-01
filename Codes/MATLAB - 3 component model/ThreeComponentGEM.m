function [t,T,alb, emm] =  ThreeComponentGEM(C_T,Q0,sigma_0,alpha_0, muFun,eps_A, eps_E, nu, tspan, y0)
%% A three component global energy budget model
% This function performs numerical simulation of the following model
% C_T T' = Q0 (1 - alpha) - sigma T^4 + mu + nu * W'
% alpha' = - eps_A * (alpha - alpha_0(T))
% sigma' = - eps_E * (simga - sigma_0(T))
% The stochastic differential equation is integrated using a standard
% forward Euler integration scheme.
%
% INPUT:
% C_T: heat capacity
% Q0: incoming solar radiation
% sigma_0: equilibrium emissivity for given temperature (function of T)
% alpha_0: equilibrium albedo for given temperature (function of T)
% muFun: the CO2 forcing (given as function of time t)
% eps_A: relaxation rate for albedo
% eps_E: relaxation rate for emissivity
% nu: variance of white noise in simulation
% tspan: specifies (approximate) times to obtain values for
% y0: initial values for all variables
%
% OUTPUT:
% t: times at which values have been saved
% T: time series for temperature T
% alb: time series for albedo alpha
% emm: time series for emissivity sigma


%% Set-up

dt = 0.001; % Time-step; 0.001 seems to work fine

%% Initial values
t = 0;
y = y0;
i = 1; % tracks amount of saved states

%% Simulation loop

while(t < tspan(end))
    
    if( t >= tspan(i) ) % save data closest to inputted set
        y_sav(i,:) = y;
        i = i + 1;
    end
    
    % Euler iteration
    y = y + dt * odefun(y,t,Q0,sigma_0,alpha_0,eps_A,eps_E,muFun,C_T,nu,dt);
    
    t = t + dt;
    
end

% It could happen that we did not save the last point when we needed to, so
% this is to prevent that from happening
if i <= length(tspan)
    y_sav(i,:) = y;
end

%% Construct the output time series
y = y_sav;

T = y(:,1);
alb = y(:,2);
emm = y(:,3);
t = tspan';



%% ODE function
% y' = f(t;...)
function dydt = odefun(y,t,Q0,sigma_0,alpha_0,eps_A,eps_E,muFun,C_T,nu,dt)
    
    % Obtain current values (simplifies notation)
    T = y(1);
    alpha = y(2);
    sigma = y(3);
    mu = muFun(t);
    
    % Deterministic part
    dydt = [1/C_T * (Q0 * (1 - alpha) + mu - sigma * T^4); ...
        -eps_A * (alpha - alpha_0(T)); ...
        -eps_E * (sigma - sigma_0(T))];
    % Add stochastic part
    dydt = dydt + nu/C_T / sqrt(dt) * randn() * [1;0;0];

