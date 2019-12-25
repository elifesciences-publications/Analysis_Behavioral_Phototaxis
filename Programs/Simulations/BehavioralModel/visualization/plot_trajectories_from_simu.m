load('/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/SimulGeorges/theta_HBO_simul.mat', 'theta')

dtheta = diff(theta,1,2);
[bout_indexes_rows, bout_indexes_columns] = find(diff(theta,1,2)~=0);

%% generate 1 population of inter-bout distances
sigmadispl = 0.6;
mudispl=0.9;
possible_lambda_values = 0:0.01:30;


%% generate corresponding (x,y) cartesian coordinates
th = theta;
Nexp = size(th,1);
Nsteps = size(th,2);
a = NaN(Nexp, Nsteps);
b = NaN(Nexp, Nsteps);
lambda = NaN(Nexp, Nsteps);
for i = 1 : Nsteps
    if i == 1
        a(:,i) = 0;
        b(:,i) = 0;
    else
        dth = th(:,i)-th(:,i-1);
        nonzero_dtheta = find(dth);
        l = zeros(Nexp,1);
        if ~isempty(nonzero_dtheta)
            l(nonzero_dtheta)...
                = pick_from_lognorm_distrib(possible_lambda_values, sigmadispl, mudispl, length(nonzero_dtheta))';
        end
        [da, db] = pol2cart(th(:,i), l);
        a(:,i) = a(:,i-1) + da;
        b(:,i) = b(:,i-1) + db;
    end
end

% finalize
virtual_roi_size = 42 ;
within_roi = (sqrt(a.^2 + b.^2) < virtual_roi_size);
th(within_roi == 0) = NaN;
meanThovertime = circ_mean(th);
Rcirc = circ_r(th);
Rproj = Rcirc.*cos(meanThovertime);

%***
figure
subplot(1,2,1)
histogram(wrapToPi(th))

subplot(1,2,2)
plot(Rproj)
title('Rproj')

% coordinates
aROI = a; 
bROI=b;
aROI(within_roi == 0)= NaN;
bROI(within_roi == 0)= NaN;

%***
figure
plot(aROI', bROI')
%%
% --- polar distribution of x ---
if ~exist('XLat','var')
    path = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/PooledData/';
    name = 'lateralized_exps.mat';
    load([path name], 'El')
    XLat = El.AngleSource;
    Xfilt = El.AngleSourceFiltered;
    FishID = El.FishN;
    medboutsperseq = median(size(XLat,2) - sum(isnan(XLat),2));
    clear El
end

Nbins=18;
colourID = 3;
fig = polar_distribution_simu_exp(XLat(:,1:medboutsperseq)+pi/2, th, Nbins, FishID, colourID);
title('escaping ROI')

fig = polar_distribution_simu_exp(XLat(:,1:medboutsperseq)+pi/2, theta, Nbins, FishID, colourID);
title('without escaping ROI')

%% save data for display
subset = 1:50;
ts = ones(subset(end),1) .* (0:0.01:(Nsteps-1)/100);

alpha = theta(subset,:);
x = aROI(subset,:);
y = bROI(subset,:);
rho = sqrt( (x).^2 + (y).^2 );

trajectories.alpha.comment = 'SIMU orientation towards source (rad)';
trajectories.alpha.values = alpha;
trajectories.rho.comment = 'distance from origin (mm)';
trajectories.rho.values = rho;
trajectories.time.comment = 'time stamp (sec)';
trajectories.time.values = ts;

save(['/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/Programs/TrajectoryDisplay/trajectories.mat'],...
    'trajectories')