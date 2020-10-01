function [ests, ests_info, ests_EV] =  perfom_estimations(DT,DR,DALB,DEMM, DALBd, DEMMd, C_T)
%% perform_estimations: gives estimates for DT_eq
% This function computes estimates for the real equilibrium warming DT_eq
% of the system. Estimates are given as function of time t, where each
% estimation value only uses data up to time t.
%
% INPUT:
% DT: time series for difference in temperature
% DR: time series for difference in radiative imbalance
% DALB: time series for difference in albedo value
% DEMIS: time series for difference ein emissivity value
% DALBd: time series for derivative of albedo
% DEMMd: time series for derivative of emisivity
% C_T: value for heat capacity (to enable computation of eigenvalues)
%
% OUTPUT:
% ests: cell containing all estimation values
% ests_info: cell containing minor clues to the method used
% ests_EV: cell containing eigenvalue estimates (if applicable)

%% Estimate 1: raw temperature time series
ests{1} = DT;
ests_info{1} = 'RAW';
ests_EV{1} = NaN;

%% Estimate 2: Linear Gregory fitting
DT_greg = NaN(1,length(DT));
EV_greg = NaN(1,length(DT));
for i = 1:length(DT)
    coeff = polyfit(DT(1:i),DR(1:i),1);
    DT_greg(i) = coeff(2) / (-coeff(1));
    EV_greg(i) = coeff(1)/C_T; % Dividy by C_T since R = C_T dT/dt
end
ests{2} = DT_greg;
ests_info{2} = 'Gregory';
ests_EV{2} = EV_greg;

%% Estimate 3: Double Gregory

% Construct a double linear function to fit the data to
f = @(x,xdata) (x(3)+x(1)*xdata) .* (xdata < - (x(4)-x(3))/(x(2)-x(1))) + ...
    (x(4)+x(2)*xdata) .* (xdata >= - (x(4)-x(3))/(x(2)-x(1)));

% To have good initial guess use linear regression on first 10 datapoints
coeff = polyfit(DT(1:10), DR(1:10),1);
lambda_1 = coeff(1);
f_1 = coeff(2);

DT_doubleGreg = NaN(1,length(DT));
EV_doubleGreg = NaN(2,length(DT));

opts1=  optimset('display','off');

for i = 1:length(DT)
    x = lsqcurvefit(f,[lambda_1, lambda_1 + 0.1, f_1, f_1], DT(1:i), DR(1:i),[],[], opts1);
    DT_doubleGreg(i) = - x(4)/x(2);
    EV_doubleGreg(:,i) = [x(1);x(2)]/C_T; % Dividy by C_T since R = C_T dT/dt
end
ests{3} = DT_doubleGreg;
ests_info{3} = 'Double Gregory';
ests_EV{3} = EV_doubleGreg;

%% Estimate 4: Linear System Fit (T & albedo)

yis = NaN(2,length(DT));
EVs = NaN(2,length(DT));

% first 5 time points do not contain enough data to fit the linear system
for i = 5:length(DT)
    Y = [DR(1:i),DALBd(1:i)];
    [n,d] = size(Y);
    X = [ones(n,1), DT(1:i), DALB(1:i)]; % column of ones for intercept
    BETA = mvregress(X,Y); % the regression
    
    % Fitted matrix A:
    A = [ BETA(2,1), BETA(3,1); ...
          BETA(2,2), BETA(3,2)];
    % Fitted vector b
    b = [ BETA(1,1); BETA(1,2) ];
    
    % equilibrium estimation
    yis(:,i) = - inv(A) * b;
    
    % computation of eigenvalues (need row operation of first row, because
    % of presence of prefactor C_T.
    A(1,:) = A(1,:)/C_T;
    lambdas = sort(eig(A));
    EVs(:,i) = lambdas;
    
end
    
ests{4} = yis(1,:);
ests_info{4} = 'Sys.Fit[T,ALB]';
ests_EV{4} = EVs;
    
%% Estimate 5: Linear System Fit (T & EMM)

yis = NaN(2,length(DT));
EVs = NaN(2,length(DT));

% first 5 time points do not contain enough data to fit the linear system
for i = 5:length(DT)
    Y = [DR(1:i),DEMMd(1:i)];
    [n,d] = size(Y);
    X = [ones(n,1), DT(1:i), DEMM(1:i)]; % column of ones for intercept
    BETA = mvregress(X,Y); % the regression
    
    % Fitted matrix A:
    A = [ BETA(2,1), BETA(3,1); ...
          BETA(2,2), BETA(3,2)];
    % Fitted vector b
    b = [ BETA(1,1); BETA(1,2) ];
    
    % equilibrium estimation
    yis(:,i) = - inv(A) * b;
    
    % computation of eigenvalues (need row operation of first row, because
    % of presence of prefactor C_T.
    A(1,:) = A(1,:)/C_T;
    lambdas = sort(eig(A));
    EVs(:,i) = lambdas;
    
end
    
ests{5} = yis(1,:);
ests_info{5} = 'Sys.Fit[T,EMM]';
ests_EV{5} = EVs;
    
%% Estimate 6: Linear System Fit (T & albedo & EMM)

yis = NaN(3,length(DT));
EVs = NaN(3,length(DT));

% first 10 time points do not contain enough data to fit the linear system
for i = 10:length(DT)
    Y = [DR(1:i),DALBd(1:i),DEMMd(1:i)];
    [n,d] = size(Y);
    X = [ones(n,1), DT(1:i), DALB(1:i),DEMM(1:i)]; % column of ones for intercept
    BETA = mvregress(X,Y); % the regression
    
    % Fitted matrix A:
    A = [ BETA(2,1), BETA(3,1), BETA(4,1); ...
          BETA(2,2), BETA(3,2), BETA(4,2); ...
          BETA(2,3), BETA(3,3), BETA(4,3)];
    % Fitted vector b
    b = [ BETA(1,1); BETA(1,2); BETA(1,3) ];
    
    % equilibrium estimation
    yis(:,i) = - inv(A) * b;
    
    % computation of eigenvalues (need row operation of first row, because
    % of presence of prefactor C_T.
    A(1,:) = A(1,:)/C_T;
    lambdas = sort(eig(A));
    EVs(:,i) = lambdas;
    
end
    
ests{6} = yis(1,:);
ests_info{6} = 'Sys.Fit[T,ALB, EMM]';
ests_EV{6} = EVs;
    







    
    
    
    
    
    
    
