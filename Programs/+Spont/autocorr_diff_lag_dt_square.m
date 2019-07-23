function [mV1, binvals] = autocorr_diff_lag_dt_square(dX,dT)
% AUTOCORRELATION of variable at different lags

variance = nanvar(dX(:));
minbin = 5;
maxbin = 7;

for dt = 1:dT
    
    Vart1 = (dX(:, 1:end-dt-1)).^2;
    Vart2 = dX(:, dt+1:end-1).^2/variance;% - dX(:, 2:end-dt); 
    
    [binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(Vart1, Vart2, minbin, maxbin);
    
    mV1 = nanmean(v2bin,2);
    stdV1 = std(v2bin,1,2);
     
end

mV1 = mV1';
stdV1 = stdV1';
binvals = binvals';

end
