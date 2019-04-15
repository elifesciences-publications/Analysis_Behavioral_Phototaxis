function[errfit, mu, sigma] = theta_measure_error_estimation(varargin)

% Estimation of Theta (X) measure error
% => gaussian fit on baseline (noise)

<<<<<<< HEAD
pathdirectory = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/Uncertainty';
=======
pathdirectory = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/Error_estimation_theta';
>>>>>>> master
filename = 'X_error_estimation';

load([pathdirectory filesep filename '.mat'], 'Err')

mu = Err.Mean;
sigma = Err.Sigma;
errfit = 0;

if nargin > 0
    x = varargin{1};
    errfit = 1/(sqrt(2*pi)*sigma)*exp( -(x - mu).^2 / (2*sigma^2));
end