function [ACnorm, sem] = xcorrMatrixRows (M)
%% biased autocorrelation on matrix rows
% input : M n*m matrix
% output :  % AC norm : normalized mean autocorrelation
            % AC 1*m vector of mean autocorrelation

% first remove starting NaNs
for i = 1 : size(M,1)       
    m = M(i,:);
    nanindex = find(diff(isnan(m))==-1, 1, 'first');
    if ~isempty(nanindex)
        if sum(isnan(m(1:nanindex))) == length(m(1:nanindex))
            m = m(nanindex+1 : end);
        end
    end
    M(i,:) = nan;
    M(i,1:length(m)) = m;
end
            
%%
% calculate autocorr on new matrix            
ac = NaN(size(M));

for i = 1:size(M,1)
    mi = M(i,:); 
    firstnan = find(isnan(mi),1, 'first');
    if isempty(firstnan)
        mo = mi;
        [xc, lag] = xcorr(mo, 'biased');
        ac(i,:) = xc( find(lag==0) : end);
    elseif firstnan > 2
        mo = mi(1:firstnan-1);
        [xc, lag] = xcorr(mo, 'biased');
        ac(i,1:firstnan-1) = xc( find(lag==0) : end);
    else
        ac(i,1:firstnan-1) = NaN;
    end
end

ac_norm = ac./ac(:,1);
ACnorm = nanmean(ac_norm,1);

sem = nanstd(ac_norm,1)./sqrt(size(ac_norm,1)-sum(isnan(ac_norm), 1));

end