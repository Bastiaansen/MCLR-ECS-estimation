close all

l1 = -2;
l2 = -0.2;

DT = 10;

a1 = -3;
a2 = - (DT+a1)

ts = 0:0.1:20;

DTs = DT + a1 * exp(l1*ts) + a2 * exp(l2*ts);
DRs = a1 * l1 * exp(l1*ts) + a2*l2*exp(l2*ts);

figure()
scatter(DTs,DRs, 'ro')
hold on


% long-time
Y = DRs;
X = DTs;


coeff = polyfit(X,Y,1);
AAA= coeff(1);
FFF = coeff(2);

plot([0, -FFF/AAA], [FFF,0], 'k-')


% short-time
Y = DRs(1:10);
X = DTs(1:10);

coeff = polyfit(X,Y,1);
AAA= coeff(1);
FFF = coeff(2);

plot([0, -FFF/AAA], [FFF,0], 'b-')

% long-time
Y = DRs(end-10:end);
X = DTs(end-10:end);


coeff = polyfit(X,Y,1);
AAA= coeff(1);
FFF = coeff(2);

plot([0, -FFF/AAA], [FFF,0], 'y-')




legend()

