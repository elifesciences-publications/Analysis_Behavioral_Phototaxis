function [ang_bout_labini, coordinatesBouti]=findAngBout(boutIndices,mang_labini, coordinates_ini)
%% Comments
% Find the angle before and after a turn bout and do the distribution per
% sequence and for the fish
%Inputs
% ind_b_a_bout: from the angular_velocity_ini function, index of the angle
% mang_labini: movmean of the angle, lab reference
%Outputs
% ang_bout_labini: bout angle, per sequence, lab reference


%% Code
[n, m] = size(boutIndices);
ang_bout_labini = nan(n, m);
coordinatesBouti = nan(n, m, 2);

for i = 1:size(boutIndices,1)
    f = find(isnan(boutIndices(i,:)),1)-1;
    if isempty(f) == 0
        ang_bout_labini(i,1:f) = mang_labini(i,boutIndices(i,1:f));
        coordinatesBouti(i,1:f,:) = coordinates_ini(boutIndices(i,1:f), :, i);
    else
        ang_bout_labini(i,:) = mang_labini(i,boutIndices(i,:));
        coordinatesBouti(i,:,:) = coordinates_ini(boutIndices(i,:), :, i);
    end
end

