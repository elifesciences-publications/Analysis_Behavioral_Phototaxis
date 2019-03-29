function [Es] = load_spont_data()

% ---
timemin = 5;
minBoutNumber = 3;
% ---

path = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/PooledData/';
name = 'spontaneous_swim.mat';

if ~exist([path name], 'file')
    [Xi, TimeBouti, xCoordi, yCoordi, FishNi] = Spont.all_experiments_pool(timemin, minBoutNumber);
    [Es] = Spont.spont_save_pooled_exps(Xi, TimeBouti, xCoordi, yCoordi, FishNi);
else
    load([path name], 'Es')
end
