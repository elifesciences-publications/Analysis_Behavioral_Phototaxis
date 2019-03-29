function[binvals, v2mean, v2std, v2eltspb] = regular_bins(Vart1, Vart2, Nbins)

binvals = min(Vart1(:)) : range(Vart1(:))/Nbins : max(Vart1(:));
v2mean = NaN(1, length(binvals)-1);
v2std = NaN(1, length(binvals)-1);
v2eltspb = NaN(1, length(binvals)-1);
for j = 2 : length(binvals)
    sel = Vart2(Vart1 > binvals(j-1) & Vart1 < binvals(j));
    v2mean(j-1) = nanmean(sel);
    v2std(j-1) = nanstd(sel);
    v2eltspb(j-1) = sum(~isnan(sel));
end
binvals = binvals(1:end-1) + mean(diff(binvals));