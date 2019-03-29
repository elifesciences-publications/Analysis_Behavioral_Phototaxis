function[A] = get_drift_for_theta(x)

% for linearly varying contrast
% iL = iR at pi/2 and -pi/2

x = wrapToPi(x);

a = 0.11;
b=0;

A = a*x + b;


