function [binvals,elts_per_bin, CbinMatrix, wturn, wfor, pturn, pfor] = probaSwitchingSides(Lum, dX, min_bins, max_bins)

dLum = diff(Lum, 1,2);

dII = dLum(:, 1:end-2)./( (Lum(:, 1:end-3) + Lum(:, 2:end-2))./2 ) ;
C = dX(:, 2:end-1).*dX(:, 3:end)./(abs(dX(:, 2:end-1)).*abs(dX(:, 3:end)));
dXresized = dX(:, 2:end);

dII = dII(:);
C = C(:);

dII(isnan(C)) = [];
C(isnan(C)) = [];

nb_elements = numel(C); % total number of elements

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

%%
% sort variable 1 ascending & sort variable 2 accordingly
% -------------------------------------------------------------------------
[~, indsort] = sort(dII);
Csorted = C(indsort);
dIIsorted = dII(indsort);
dXsorted = dXresized(indsort);
bins = [0 ; cumsum(ones(nbins, 1)*elts_per_bin)]+1;

CbinMatrix = NaN(nbins-1, elts_per_bin);
binvals = NaN(nbins-1,1);
wturn = NaN(nbins-1,1);
wfor = NaN(nbins-1,1);
pturn = NaN(nbins-1,1);
pfor = NaN(nbins-1,1);
for j = 2 : nbins
    CbinMatrix(j-1, :) = Csorted(bins(j-1):bins(j)-1)';
    dXbin = dXsorted(bins(j-1):bins(j)-1)';
    [wturn(j-1), wfor(j-1), pturn(j-1), pfor(j-1)] = Gauss2custom_simple(dXbin);
    binvals(j-1) = mean(dIIsorted(bins(j-1):bins(j)-1));
end

q = 0.5 - pi/4*mean(CbinMatrix, 2).*( 1+ (pfor./pturn.*wfor./wturn).^2 );

%***
figure
plot(binvals, q);

%% stuff to display

disp(['nb elements : ' num2str(nb_elements)])
disp(['number of bins ' num2str(nbins-1)])
disp(['nb of elements per bin ' num2str(elts_per_bin)])
disp(['omitted elements ' num2str(minmod)])
