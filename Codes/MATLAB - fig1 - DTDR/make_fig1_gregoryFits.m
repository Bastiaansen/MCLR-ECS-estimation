% This script uses the standard 'Gregory' linear regression on DT and DR to
% obtain different estimates, based on the years used in the regression

%% Linear Regression years 1-i_end
coeff = polyfit(DT(1:i_end),DR(1:i_end),1);
x_intercept = coeff(2) / (-coeff(1));
y_intercept = coeff(2);

plot([0 2*x_intercept], [1,-1] * y_intercept, 'b-', 'linewidth', 2.0)

scatter(x_intercept, 0, 'b*')

%% Linear Regression years 20-i_end
coeff = polyfit(DT(20:i_end),DR(20:i_end),1);
x_intercept = coeff(2) / (-coeff(1));
y_intercept = coeff(2);

plot([0 2*x_intercept], [1,-1] * y_intercept, 'y-', 'linewidth', 2.0)

scatter(x_intercept, 0, 'y*')

%% Linear Regression years 1-20
coeff = polyfit(DT(1:20),DR(1:20),1);
x_intercept = coeff(2) / (-coeff(1));
y_intercept = coeff(2);

plot([0 2*x_intercept], [1,-1] * y_intercept, 'k-', 'linewidth', 2.0)

scatter(x_intercept, 0, 'k*')

%% Linear Regression last 100 years

coeff = polyfit(DT(end-100:end),DR(end-100:end),1);
x_intercept = coeff(2) / (-coeff(1));
y_intercept = coeff(2);

plot([0 2*x_intercept], [1,-1] * y_intercept, 'k--', 'linewidth', 2.0)

scatter(x_intercept, 0, 'k*')