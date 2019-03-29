function [lumBout, angleSourceBout, angleLabSmBout, coordinatesBout] = findAngleLumBout(boutIndices, luminosity, angleSourceSm, angleLabSm, coordinates) 

[n, m] = size(boutIndices);
lumBout = nan(n, m);
angleSourceBout = nan(n, m);
coordinatesBout = nan(n, m, 2);
angleLabSmBout = nan(n, m);

for i = 1:size(boutIndices,1)
    % if 0 in bout indices 
    % f = find(ind_b_a_bout(i,:)==0, 1) - 1; 
    % if NaNs
    f = find(isnan(boutIndices(i,:)), 1) - 1;
    if ~isempty(f)
        coordinatesBout(i,1:f,:) = coordinates( boutIndices(i,1:f), :, i);
        lumBout(i,1:f) = luminosity(i,boutIndices(i,1:f));
        angleSourceBout(i,1:f) = angleSourceSm(i,boutIndices(i,1:f));
        angleLabSmBout(i,1:f) = angleLabSm(i,boutIndices(i,1:f));
    else
        coordinatesBout(i,:,:) = coordinates(boutIndices(i,:),:,i);
        lumBout(i,:) = luminosity(i,boutIndices(i,:));
        angleSourceBout(i,:) = angleSourceSm(i,boutIndices(i,:));
        angleLabSmBout(i,:) = angleLabSm(i,boutIndices(i,:));
    end
end