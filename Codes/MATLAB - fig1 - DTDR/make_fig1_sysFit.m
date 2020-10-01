% This script uses the system fit linear regression with DY = [DR, DALBd,
% DEMMd] and DX = [DT, DALB, DEMM] to infer the equilibrium warming

%% System Fit

DY = [DR(1:i_end), ALBd(1:i_end), EMMd(1:i_end)];

[n,d] = size(DY);
% Artifact of matlab implementation: it cannot deal well with very small
% values, so input has been scaled to be O(1); since it is a linear system,
% and no additional regularization is applied this does not skew the
% results.
DX = [ones(n,1), DT(1:i_end), 10^2*DALB(1:i_end), 10^9*DEMM(1:i_end)];

BETA = mvregress(DX,DY);

A = [ BETA(2,1), BETA(3,1), BETA(4,1); ...
        BETA(2,2), BETA(3,2), BETA(4,2); ...
        BETA(2,3), BETA(3,3), BETA(4,3)];
b = [ BETA(1,1); BETA(1,2);BETA(1,3) ];
yis = - inv(A) * b;

%% PLOTTING


% MUCH SMOOTHING
% (that is needed here, since we need 3D input, but we can only infer that
% from the initial data; so without smoothing of that input, the regression
% seems to be overfitting, but that is only because of the data. In
% reality, only one path through the fitted linear 3D object will be taken,
% but there is no way to reconstruct that without using the input data --
% in contrast to the other methods)
k = 50;

% DTs = movmean(DT,k);
% DALBs = movmean(DALB,k)*10^2;
% DEMMs = movmean(DEMM,k)*10^9;

k1 = 10;
k2 = 21;
k3 = 100;
k4 = 500;

DTs = [movmean(DT(1:10),k1); movmean(DT(11:100),k2); movmean(DT(101:1000), k3); movmean(DT(1001:end),k4)];
DALBs = 10^2*[movmean(DALB(1:10),k1); movmean(DALB(11:100),k2); movmean(DALB(101:1000),k3); movmean(DALB(1001:end),k4)];
DEMMs = 10^9*[movmean(DEMM(1:10),k1); movmean(DEMM(11:100),k2); movmean(DEMM(101:1000),k3); movmean(DEMM(1001:end),k4)];

for i = 1:length(DTs)
    DXi = [DTs(i);DALBs(i);DEMMs(i)];
    DRi = A * DXi + b;
    DRs(i) = DRi(1);
end

saddleBrown = [139,	69,	19]/255*1.1;	
plot(DTs,DRs, 'color', 'm', 'linewidth', 3.0)
scatter(yis(1), 0, '*', 'MarkerEdgeColor', 'm')