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

%% Trajectories
xCoord0 = (xCoord-xCoord(:,1))/pxmm;
yCoord0 = (yCoord-yCoord(:,1))/pxmm;
rho = sqrt( (xCoord0).^2 + (yCoord0).^2 );

x = xCoord0.*cos(-wrapToPi(XLat(:,1))+Xlab(:,1)) - yCoord0.*sin(-wrapToPi(XLat(:,1)) + Xlab(:,1));
y = xCoord0.*sin(-wrapToPi(XLat(:,1))+Xlab(:,1)) + yCoord0.*cos(-wrapToPi(XLat(:,1)) + Xlab(:,1));
plot(x',y','.-')
title('source to the top')