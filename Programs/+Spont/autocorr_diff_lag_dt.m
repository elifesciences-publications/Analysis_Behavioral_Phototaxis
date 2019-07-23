function [fig, aclag1] = autocorr_diff_lag_dt(dX,dT)
AUTOCORRELATION of variable at different lags

minbin = 7;
maxbin = 9;

***
fig = figure;
hold on

for dt = 1:dT
    
    Vart1 = (dX(:, 1:end-dt-1));
    Vart2 = dX(:, dt+1:end-1);% - dX(:, 2:end-dt); 
    
    [binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(Vart1, Vart2, minbin, maxbin);
    
    mV1 = nanmean(v2bin,2);
    stdV1 = std(v2bin,1,2);
    
    %***
    hold on
    c = dt/(dT+0.1);
    errorbar(binvals(1:end-1), mV1(1:end-1), stdV1(1:end-1)/sqrt(elts_per_bin),...
         '-', 'DisplayName', ['p = ' num2str(dt)], 'Linewidth', 1.5, 'Color', [0 c c*0.8])
     drawnow
     pause(0.5)
     
end
xlim([0 pi/6])
ax = gca;
ax.XScale = 'lin';
xticks auto
xticklabels auto

xlabel('<|\delta\theta_{n-1}|>_{bias corrected}')
ylabel('dunno')
legend

ax.FontSize = 14;

end
