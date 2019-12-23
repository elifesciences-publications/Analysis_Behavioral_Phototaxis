function [ACnorm, sem, ac_fish] = xcorrMatrixRowsPerFish (M, FishID)
%% biased autocorrelation on matrix rows
% input : M n*m matrix
% output :  % AC norm : normalized mean autocorrelation
            % AC 1*m vector of mean autocorrelation

fish = unique(FishID);
ac_fish = NaN(length(fish), size(M,2));

for i = 1:length(fish)
    seqafish = (FishID == i);
    fish = M(seqafish,:);
    only_nans_behind = find(diff(sum(isnan(fish),1)), 1, 'last');
    if sum(isnan(fish(:,only_nans_behind+1)),1) == sum(seqafish)
        fish = fish(:,1:only_nans_behind);
    end
    [acfish, ~] = xcorrMatrixRows (fish);
    len = length(acfish);
    ac_fish(i,1:len) = acfish;
end

ACnorm = nanmean(ac_fish,1);

sem = nanstd(ac_fish,1)./sqrt(size(ac_fish,1)-sum(isnan(ac_fish), 1));

end