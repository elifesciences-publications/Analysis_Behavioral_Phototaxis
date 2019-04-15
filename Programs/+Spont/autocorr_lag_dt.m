function [fig, aclag1] = autocorr_lag_dt(dXissbiais,dT)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%***
fig = figure;
for dt = 1:dT
    
    Vart1 = (dXissbiais(:, 1:end-dt));
    Vart2 = (dXissbiais(:, dt+1:end)); 
    
    minbin = 18;
    maxbin = 30;
    [binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(Vart1, Vart2, minbin, maxbin);
    
    mV1 = nanmean(v2bin,2);
    stdV1 = std(v2bin,1,2);
    
    %***
    hold on
    errorbar(binvals(1:end-1), mV1(1:end-1), stdV1(1:end-1)/sqrt(elts_per_bin),...
         '-', 'DisplayName', ['AC lag ' num2str(dt)], 'Linewidth', 1.5, 'Color', [dt/5 dt/5 dt/5])
end
xticks([ -pi/3, 0, pi/3])
xticklabels({'-\pi/3', '0', '\pi/3'})
yticks([-0.2:0.1:0.2])
xlim([-pi/2 pi/2])
xlabel('<\delta\theta_n_-_1>_{bias corrected}')
ylabel('<\delta\theta_n>_{bias corrected}')
ax = gca;
ax.FontSize = 14;

aclag1 = (Vart1(:).*Vart2(:)./abs(Vart1(:).*Vart2(:)));
aclag1 = nanmean(aclag1)

end

