function [binvals, elts_per_bin, v2bin, sorted_indices, binedges] = BinsWithEqualNbofElements(v1, v2, min_bins, max_bins)

v1 = v1(:);
v2 = v2(:);

v1(isnan(v2)) = [];
v2(isnan(v2)) = [];

v2sq = v2.^2;

nb_elements = numel(v2); % total number of elements

%%
% find number of bins so that each bin has the same number of values
% -------------------------------------------------------------------------
nbins = [];
minmod = inf;
for j = min_bins:max_bins
    if mod(nb_elements, j) < minmod
        nbins = j;
        minmod = mod(nb_elements, j);
    end
end
if isempty(nbins)
    warning('check bins')
    return
end
elts_per_bin = (nb_elements - minmod)/nbins;

disp(['nb elements : ' num2str(nb_elements)])
disp(['number of bins ' num2str(nbins)])
disp(['nb of elements per bin ' num2str(elts_per_bin)])
disp(['omitted elements ' num2str(minmod)])

%%
% sort variable 1 ascending & sort variable 2 accordingly
% -------------------------------------------------------------------------
nbins = nbins+1;
[~, indsort] = sort(v1);
v2sorted = v2(indsort);
v2sqsorted = v2sq(indsort);
v1sorted = v1(indsort);
bins = [0 ; cumsum(ones(nbins, 1)*elts_per_bin)]+1;
v2bin = NaN(nbins-1, elts_per_bin);
sorted_indices = NaN(nbins-1, elts_per_bin);
v2sqbinMatrix = NaN(nbins-1, elts_per_bin);
binvals = NaN(nbins-1,1);
binedges = NaN(nbins, 1); 
binedges(1) = v1sorted(1);
for j = 2 : nbins
    v2bin(j-1, :) = v2sorted(bins(j-1):bins(j)-1)';
    sorted_indices(j-1, :) = indsort(bins(j-1):bins(j)-1)';
    v2sqbinMatrix(j-1, :) = v2sqsorted(bins(j-1):bins(j)-1)';
    binvals(j-1) = nanmean(v1sorted(bins(j-1):bins(j)-1));
    binedges(j) = v1sorted(bins(j)-1);
end
v1sorted(isnan(v1sorted))=[];
binedges(end) = v1sorted(end);
