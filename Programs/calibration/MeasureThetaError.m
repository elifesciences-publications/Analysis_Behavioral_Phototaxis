% estimate theta measurement error
timemin = 3;
EXP=1;
% choose experiment
if EXP == 1
    [Dates, Fish, exptype, ~, intmax] = Temp.chooseExpType('exp60');
elseif EXP == 2
    [Dates, Fish, exptype, ~, intmax] = Temp.chooseExpType('exp30');
elseif EXP ==3
    [Dates, Fish, exptype, ~, intmax] = Temp.chooseExpType('sin60');
end

A = [];
Mu = [];
Sigma = [];
SE = [];
count = 1;
for p = 1: length(Dates)
    for q = 1: length(Fish(p,:))
        if iscellstr(Fish(p,q)) == 0
            continue
        end
        
        fish = char(Fish(p,q));
        date = char(Dates(p));
        
        load(['/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Data/'...
            exptype '/' intmax '/' date '/' fish '/data' fish '.mat']);
        
        % --- Remove the short sequences ---
        [angleLab_ini, coordinates_ini, framerateInfo_ini] = remove_sequence_ini(D, timemin);
        framerate = framerateInfo_ini(:,3);
        
        % --- Get the cumulative angle ---
        %     (lab reference, the angle is initially wrapped to 360) 
        [angleLabCum_ini] = angle_cum(angleLab_ini, framerateInfo_ini);
        
        % --- Correct for potential aberrant bouts (head/tail misplacement)
        angleLabCumCorr_ini = find_aberrant_bouts(angleLabCum_ini, coordinates_ini, 0);
        
        % --- Spot real bouts
        fig = 0;
        boutIndices_ini = BoutSpot(coordinates_ini, angleLabCumCorr_ini, framerateInfo_ini, fig);    
        
        angle_for_error = angleLabCumCorr_ini;
        
        fig = 1;
        [a, mu, sigma, se] = angle_measurement_error(angle_for_error, boutIndices_ini, framerate, fig);
        A(count) = a;
        Mu(count) = mu;
        Sigma(count) = sigma;
        SE(count) = se;
        count = count+1;
    end
end

mu = mean(Mu);
sigma = mean(Sigma);
a = mean(A);

x = [-0.2:0.01:0.2];
errfit = 1/(sqrt(2*pi)*sigma)*exp( -(x - mu).^2 / (2*sigma^2));
plot(x, errfitnorm, 'k', 'LineWidth', 2)

%--- store this sh.. ---
pathdirectory = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/Error_estimation';
filename = 'X_error_estimation';

Err.comment = 'Error estimation from fits on trajectories from N fish';
Err.N = count-1;
Err.Mean = mu;
Err.Sigma = sigma;
Err.individual.comment = 'Error estimation from fit on different fish trajectories';
Err.individual.mean = Mu;
Err.individual.sigma = Sigma;
Err.individual.standdev = SE;

save([pathdirectory filesep filename '.mat'], 'Err')


