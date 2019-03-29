% combine temporal phototaxis data

%% === Load data ===
Temp.Load

boutsperseq = size(X,2)-sum(isnan(X),2);
mednbouts = median(boutsperseq);
%% Plot sequence biases for individual fish

% ***
figure
plot( Bias, nanmean(dX,2), 'o')
xlabel('bias per sequence')
ylabel('mean bias per fish')
%% fit IBI
fitIBI = fitdist(IBI(:), 'Gamma');
IBI(IBI> fitIBI.b*3) = NaN;

%% <dX> vs X or I
Temp.OrientationBias_dXvsX

%% Variance 
Temp.Variance

%% Inter-bout interval
Temp.InterBoutInterval

%% Distances and velocities between bouts vs I & dI/I
Temp.DistanceVelocity

%% Autocorrelation
Temp.Autocorrelation

%% Distribution analysis : double gaussian fits
Temp.dXDistribution    