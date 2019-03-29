function[AC] = dXautocorrelation(dX, fig)
% dXautocorrelation
% autocorrelations of dX

% --- calculate autocorrelation ---
ac = NaN(size(dX));
for i = 1:size(dX,1)
    mi = dX(i,:);
    firstnan = find(isnan(mi),1, 'first');
    if isempty(firstnan)
        mo = mi;
        xc = xcorr(mo, 'biased'); 
        ac(i,:) = xc(end-(length(xc)+1)/2:end);
    elseif firstnan > 2
        mo = mi(1:firstnan-1);
        xc = xcorr(mo, 'biased');
        ac(i,1:firstnan-1) = xc(firstnan-1:end); % fin-2
    else
        ac(i,1:firstnan-1) = NaN;
    end
end
AC = nanmean(ac,1);

if fig
    % ***
    figure
    plot(AC/max(AC))
    xlim([1 50])
    grid on
end