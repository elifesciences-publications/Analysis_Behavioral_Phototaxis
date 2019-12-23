%% Extract autocorrelative parameters on bouts from Lateralized (stereo-phototaxis) data
% for lag = 1
% and ACF (lag = 0 : inf)

%% ::: Autocorrelation in reinforcement/inhibition conditions :::
dX = dXLssbiais;
C = DIlr;

binsdXdX = 9;
binspflip = 7;
Lat.auto_co_reinforcement(dX, dXissbiais, C, binsdXdX, binspflip, 0, 'absolute')

%%